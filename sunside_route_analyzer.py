# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# SUNSIDE ROUTE ANALYZER
# REPOSITORY: https://github.com/cassianorcarneiro/sunside-route-analyzer
# CASSIANO RIBEIRO CARNEIRO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

"""Estimate which side of a vehicle receives less direct sunlight along a driving route.

Pipeline:
    1. Resolve places (city names via Nominatim, or raw lat/lon).
    2. Download the drivable OSM graph inside a corridor around the waypoints.
    3. Route leg by leg (origin -> via_1 -> ... -> destination).
    4. Split the route into its native OSM polyline segments (no resampling).
    5. Assign each segment a time from OSM edge travel times, rescaled to the user's arrival time.
    6. Compute solar position and clear-sky DNI per segment (vectorized pvlib).
    7. Integrate direct-beam radiant exposure on each vertical side window over time.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import pvlib
from geopy.exc import GeopyError
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim
from pyproj import Geod
from shapely.geometry import LineString

logger = logging.getLogger("sunside_route_analyzer")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Data types
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class LatLon(NamedTuple):
    """A WGS84 coordinate in decimal degrees."""

    lat: float
    lon: float


LatLonLike = tuple[float, float] | LatLon
Place = str | LatLonLike          # "Juiz de Fora, MG, Brasil", "-21.76,-43.35" or (lat, lon)
TimeLike = str | pd.Timestamp

_LATLON_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)\s*$")


@dataclass(frozen=True)
class ResolvedPlace:
    """A place after geocoding: the coordinate used for routing plus a human-readable label."""

    query: str
    label: str
    coord: LatLon


@dataclass(frozen=True)
class _Polyline:
    """Route as arrays of native OSM segments (one entry per consecutive vertex pair)."""

    lon1: np.ndarray
    lat1: np.ndarray
    lon2: np.ndarray
    lat2: np.ndarray
    length_m: np.ndarray
    heading_deg: np.ndarray
    duration_s: np.ndarray


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Core analyzer
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class SunSideRouteAnalyzer:
    """Integrates direct-beam solar exposure on each side window of a vehicle along an OSM route.

    Physical model:
        Each side window is a vertical plane whose outward normal is perpendicular to the heading.
        cos(incidence) on the left window = cos(el) * sin(heading - azimuth).
        Exposure [Wh/m^2] = sum over segments of DNI * max(0, cos(incidence)) * dt / 3600.
        The right window has the opposite sign.
    """

    _MIN_ENDPOINT_DISTANCE_M = 1000.0
    _SNAP_WARN_DISTANCE_M = 2000.0
    _TIMING_SCALE_WARN = (0.6, 1.8)
    # Free-flow speeds for untagged edges (km/h). maxspeed coverage is sparse in much of Brazil.
    _HWY_SPEEDS_KPH = {
        "motorway": 100, "motorway_link": 60, "trunk": 90, "trunk_link": 50,
        "primary": 80, "primary_link": 40, "secondary": 60, "secondary_link": 40,
        "tertiary": 50, "tertiary_link": 30, "unclassified": 40, "residential": 30,
        "living_street": 10, "service": 20,
    }
    _FALLBACK_SPEED_KPH = 40.0
    _SETTLEMENT_TYPES = frozenset({"city", "town", "village", "municipality", "hamlet"})
    _DEFAULT_GEOCODER_ZOOM_LEVELS = (10, 12, 8)
    _MUNICIPAL_KEYS = (
        "city", "town", "village", "municipality",
        "city_district", "district", "borough",
        "county", "state_district",
    )

    def __init__(
        self,
        tz: str = "America/Sao_Paulo",
        network_type: str = "drive",
        weight: str = "travel_time",
        corridor_buffer_m: float = 8000.0,
        weighting: str = "clearsky",
        balance_threshold: float = 0.15,
        use_osmnx_cache: bool = True,
        log_console: bool = False,
        geocoder_user_agent: str = "sunside-route-analyzer",
        geocoder_language: str = "pt-BR",
        geocoder_min_delay_s: float = 1.1,
    ) -> None:
        """
        Args:
            tz:                 IANA timezone used to localize naive timestamps.
            network_type:       OSMnx network type.
            weight:             Edge attribute minimized by the router ("travel_time" or "length").
            corridor_buffer_m:  Buffer distance on EACH side of the waypoint polyline.
            weighting:          "clearsky" weights incidence by Ineichen clear-sky DNI (W/m^2).
                                "geometric" uses DNI = 1 (pure geometry, unitless).
            balance_threshold:  |balance| below this returns "either" instead of a side.
        """
        if corridor_buffer_m <= 0:
            raise ValueError("corridor_buffer_m must be positive.")
        if weighting not in ("clearsky", "geometric"):
            raise ValueError("weighting must be 'clearsky' or 'geometric'.")
        if not 0.0 <= balance_threshold < 1.0:
            raise ValueError("balance_threshold must be in [0, 1).")

        self.tz = tz
        self.network_type = network_type
        self.weight = weight
        self.corridor_buffer_m = float(corridor_buffer_m)
        self.weighting = weighting
        self.balance_threshold = float(balance_threshold)
        self.geod = Geod(ellps="WGS84")

        ox.settings.use_cache = bool(use_osmnx_cache)
        ox.settings.log_console = bool(log_console)

        self.geocoder_language = geocoder_language
        self._geolocator = Nominatim(user_agent=geocoder_user_agent, timeout=10)
        self._geocode = RateLimiter(self._geolocator.geocode, min_delay_seconds=geocoder_min_delay_s)
        self._reverse = RateLimiter(self._geolocator.reverse, min_delay_seconds=geocoder_min_delay_s)
        self._geocode_cache: dict[str, ResolvedPlace] = {}

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def resolve_place(self, place: Place, name: str = "place") -> ResolvedPlace:
        """Turn a city name, a 'lat,lon' string or a (lat, lon) pair into a ResolvedPlace.

        Names are geocoded with Nominatim. Settlement results (city/town/village) are
        preferred over streets or POIs with the same name.
        """
        if isinstance(place, str):
            m = _LATLON_RE.match(place)
            if m:
                coord = self._coerce_latlon((float(m.group(1)), float(m.group(2))), name)
                return ResolvedPlace(place, f"{coord.lat:.5f},{coord.lon:.5f}", coord)
            return self._geocode_name(place, name)

        coord = self._coerce_latlon(place, name)
        return ResolvedPlace(str(tuple(coord)), f"{coord.lat:.5f},{coord.lon:.5f}", coord)

    def analyze(
        self,
        origin: Place,
        destination: Place,
        depart_time: TimeLike,
        arrive_time: TimeLike | None = None,
        vias: list[Place] | None = None,
        save_html_map: bool = True,
        html_path: str | Path = "route.html",
        zoom_start: int = 9,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        """Run the full pipeline.

        Args:
            origin, destination: City name ("Bicas, MG, Brasil"), "lat,lon" string or (lat, lon).
            depart_time:         Departure time.
            arrive_time:         Optional arrival time. If omitted, OSM travel times are used as-is.
                                 If given, OSM travel times are rescaled uniformly to match it.
            vias:                Intermediate waypoints. The route passes through each one in order.

        Returns:
            (segments_dataframe, summary_dict)
        """
        places = [self.resolve_place(origin, "origin")]
        places += [self.resolve_place(v, f"via[{i}]") for i, v in enumerate(vias or [])]
        places.append(self.resolve_place(destination, "destination"))
        waypoints = [p.coord for p in places]
        for p in places:
            logger.info("Resolved %r -> %s (%.5f, %.5f)", p.query, p.label, p.coord.lat, p.coord.lon)

        self._check_endpoint_distance(waypoints[0], waypoints[-1])
        depart, arrive = self._normalize_times(depart_time, arrive_time)

        logger.info("Downloading corridor graph (buffer=%.0f m per side)...", self.corridor_buffer_m)
        G = self._build_corridor_graph(waypoints)
        # Only edges WITHOUT a maxspeed tag use hwy_speeds; tagged edges keep their own value.
        G = ox.add_edge_speeds(G, hwy_speeds=self._HWY_SPEEDS_KPH, fallback=self._FALLBACK_SPEED_KPH)
        G = ox.add_edge_travel_times(G)

        route = self._route_through_waypoints(G, waypoints)
        poly = self._route_polyline(G, route)

        routed_s = float(poly.duration_s.sum())
        if arrive is None:
            scale = 1.0
            arrive = (depart + pd.to_timedelta(routed_s, unit="s")).round("s")
        else:
            scale = (arrive - depart).total_seconds() / routed_s
            lo, hi = self._TIMING_SCALE_WARN
            if not lo <= scale <= hi:
                logger.warning(
                    "Your trip duration is %.2fx the OSM estimate (%.0f min). "
                    "Check the times or the route.", scale, routed_s / 60.0,
                )

        df = self._compute_segments(poly, depart, scale)
        summary = self._summarize(df, places, depart, arrive, routed_s, scale)

        if save_html_map:
            summary["route_map_html"] = self.save_route_map_html(
                df, places, out_html=html_path, zoom_start=zoom_start
            )

        logger.info("Done: %d segments, %.1f km.", len(df), summary["total_distance_m"] / 1000.0)
        return df, summary

    def municipality_from_geopy(self, coordinate: LatLonLike, language: str | None = None) -> str | None:
        """Reverse-geocode a municipality name (network request). Returns None on failure."""
        coord = self._coerce_latlon(coordinate, "coordinate")
        for zoom in self._DEFAULT_GEOCODER_ZOOM_LEVELS:
            try:
                loc = self._reverse(
                    (coord.lat, coord.lon), exactly_one=True,
                    language=language or self.geocoder_language, zoom=zoom, addressdetails=True,
                )
            except GeopyError as exc:
                logger.debug("Reverse geocoding failed at zoom=%d: %s", zoom, exc)
                continue
            if not loc:
                continue
            addr = loc.raw.get("address", {})
            for key in self._MUNICIPAL_KEYS:
                if addr.get(key):
                    return addr[key]
        return None

    def save_route_map_html(
        self,
        df: pd.DataFrame,
        places: list[ResolvedPlace],
        out_html: str | Path = "route.html",
        zoom_start: int = 9,
    ) -> str:
        """Folium map with the route colored by the recommended (shadier) side per stretch."""
        colors = {"left": "#1f77b4", "right": "#ff7f0e", "none": "#7f7f7f"}
        labels = {"left": "sit left", "right": "sit right", "none": "no direct sun"}

        lats = [p.coord.lat for p in places]
        lons = [p.coord.lon for p in places]
        m = folium.Map(location=[float(np.mean(lats)), float(np.mean(lons))], zoom_start=zoom_start)

        # Group consecutive segments with the same side into one polyline.
        run_id = (df["side_less_sun"] != df["side_less_sun"].shift()).cumsum()
        for _, run in df.groupby(run_id, sort=False):
            side = run["side_less_sun"].iloc[0]
            pts = list(zip(run["lat1"], run["lon1"])) + [(run["lat2"].iloc[-1], run["lon2"].iloc[-1])]
            folium.PolyLine(
                pts, color=colors[side], weight=5, opacity=0.85,
                tooltip=f"{labels[side]} | {run['time'].iloc[0]:%H:%M}",
            ).add_to(m)

        folium.Marker([lats[0], lons[0]], tooltip=f"Origin: {places[0].label}",
                      icon=folium.Icon(color="green")).add_to(m)
        folium.Marker([lats[-1], lons[-1]], tooltip=f"Destination: {places[-1].label}",
                      icon=folium.Icon(color="red")).add_to(m)
        for i, p in enumerate(places[1:-1], start=1):
            folium.CircleMarker([p.coord.lat, p.coord.lon], radius=6, fill=True,
                                tooltip=f"Via {i}: {p.label}").add_to(m)

        legend = (
            '<div style="position:fixed;bottom:24px;left:24px;z-index:9999;background:white;'
            'padding:8px 12px;border-radius:6px;font:13px sans-serif;box-shadow:0 1px 4px #0004">'
            + "<br>".join(f'<span style="color:{colors[k]}">&#9632;</span> {labels[k]}' for k in colors)
            + "</div>"
        )
        m.get_root().html.add_child(folium.Element(legend))

        out_html = str(out_html)
        m.save(out_html)
        return out_html

    # -------------------------------------------------------------------------
    # Geocoding, validation, time handling
    # -------------------------------------------------------------------------

    def _geocode_name(self, query: str, name: str) -> ResolvedPlace:
        key = query.strip().lower()
        if key in self._geocode_cache:
            return self._geocode_cache[key]

        try:
            results = self._geocode(
                query, exactly_one=False, limit=5,
                addressdetails=True, language=self.geocoder_language,
            )
        except GeopyError as exc:
            raise RuntimeError(f"Geocoding failed for {name} {query!r}: {exc}") from exc
        if not results:
            raise ValueError(f"No geocoding result for {name} {query!r}. Try 'City, State, Country'.")

        best = next(
            (r for r in results
             if r.raw.get("addresstype") in self._SETTLEMENT_TYPES
             or (r.raw.get("class") == "place" and r.raw.get("type") in self._SETTLEMENT_TYPES)),
            None,
        )
        if best is None:
            best = results[0]
            logger.warning("%s %r did not match a city/town; using %r.", name, query, best.address)

        resolved = ResolvedPlace(
            query=query,
            label=best.address,
            coord=self._coerce_latlon((best.latitude, best.longitude), name),
        )
        self._geocode_cache[key] = resolved
        return resolved

    @staticmethod
    def _coerce_latlon(value: LatLonLike | None, name: str) -> LatLon:
        if value is None:
            raise ValueError(f"{name} must be provided.")
        try:
            lat, lon = value  # works for LatLon and plain tuples
            lat_f, lon_f = float(lat), float(lon)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be a numeric (lat, lon) pair.") from exc
        if not (-90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0):
            raise ValueError(f"{name} is out of range: ({lat_f}, {lon_f}).")
        return LatLon(lat_f, lon_f)

    def _check_endpoint_distance(self, a: LatLon, b: LatLon) -> None:
        _, _, dist_m = self.geod.inv(a.lon, a.lat, b.lon, b.lat)
        if dist_m < self._MIN_ENDPOINT_DISTANCE_M:
            raise ValueError(
                f"Origin and destination are only {dist_m:.0f} m apart "
                f"(minimum {self._MIN_ENDPOINT_DISTANCE_M:.0f} m)."
            )

    def _to_tz(self, t: TimeLike, name: str) -> pd.Timestamp:
        try:
            ts = pd.Timestamp(t)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid datetime for {name}: {t!r}") from exc
        if pd.isna(ts):
            raise ValueError(f"Invalid datetime for {name}: {t!r}")
        return ts.tz_localize(self.tz) if ts.tzinfo is None else ts.tz_convert(self.tz)

    def _normalize_times(
        self, depart_time: TimeLike, arrive_time: TimeLike | None
    ) -> tuple[pd.Timestamp, pd.Timestamp | None]:
        depart = self._to_tz(depart_time, "depart_time")
        if arrive_time in (None, ""):
            return depart, None
        arrive = self._to_tz(arrive_time, "arrive_time")
        if arrive <= depart:
            raise ValueError("arrive_time must be later than depart_time.")
        return depart, arrive

    # -------------------------------------------------------------------------
    # Graph and routing
    # -------------------------------------------------------------------------

    def _build_corridor_graph(self, waypoints: list[LatLon]) -> nx.MultiDiGraph:
        """Download the OSM graph inside a buffer around the waypoint polyline.

        The buffer is computed in the local UTM zone. Web Mercator (EPSG:3857) inflates
        distances by 1/cos(lat), so a "8 km" buffer there is ~7.4 km true at 21.5°S.
        """
        line = gpd.GeoSeries([LineString([(w.lon, w.lat) for w in waypoints])], crs="EPSG:4326")
        utm = line.estimate_utm_crs()
        corridor = line.to_crs(utm).buffer(self.corridor_buffer_m).to_crs("EPSG:4326").iloc[0]
        return ox.graph_from_polygon(corridor, network_type=self.network_type)

    def _nearest_node(self, G: nx.MultiDiGraph, p: LatLon) -> tuple[int, float]:
        """Haversine nearest node on an unprojected graph. Returns (node_id, distance_m)."""
        ids = np.fromiter(G.nodes, dtype=np.int64, count=G.number_of_nodes())
        lon = np.radians([G.nodes[n]["x"] for n in ids])
        lat = np.radians([G.nodes[n]["y"] for n in ids])
        plat, plon = np.radians(p.lat), np.radians(p.lon)
        a = np.sin((lat - plat) / 2) ** 2 + np.cos(plat) * np.cos(lat) * np.sin((lon - plon) / 2) ** 2
        d = 2 * 6_371_008.8 * np.arcsin(np.sqrt(a))
        k = int(np.argmin(d))
        return int(ids[k]), float(d[k])

    def _route_through_waypoints(self, G: nx.MultiDiGraph, waypoints: list[LatLon]) -> list[int]:
        """Shortest path through every waypoint in order (one leg per consecutive pair)."""
        nodes = []
        for i, wp in enumerate(waypoints):
            n, d = self._nearest_node(G, wp)
            if d > self._SNAP_WARN_DISTANCE_M:
                logger.warning("Waypoint %d snapped to a node %.1f km away.", i, d / 1000.0)
            nodes.append(n)

        route = [nodes[0]]
        for leg, (a, b) in enumerate(itertools.pairwise(nodes)):
            if a == b:
                continue
            try:
                path = nx.shortest_path(G, a, b, weight=self.weight)
            except nx.NetworkXNoPath as exc:
                raise RuntimeError(
                    f"No drivable path for leg {leg}. Increase corridor_buffer_m or add vias."
                ) from exc
            route.extend(path[1:])
        if len(route) < 2:
            raise RuntimeError("Route collapsed to a single node.")
        return route

    def _route_polyline(self, G: nx.MultiDiGraph, route: list[int]) -> _Polyline:
        """Expand the node route into native OSM vertex-to-vertex segments.

        Each edge geometry is oriented u -> v explicitly. A reversed geometry would flip
        the heading by 180° and swap left/right on that edge.
        Edge travel time is split over its sub-segments proportionally to length.
        """
        lon1, lat1, lon2, lat2, edge_idx, edge_tt = [], [], [], [], [], []
        for k, (u, v) in enumerate(itertools.pairwise(route)):
            data = min(G[u][v].values(), key=lambda d: d.get(self.weight, np.inf))
            ux, uy = G.nodes[u]["x"], G.nodes[u]["y"]
            vx, vy = G.nodes[v]["x"], G.nodes[v]["y"]
            if "geometry" in data:
                coords = list(data["geometry"].coords)
                d_first = (coords[0][0] - ux) ** 2 + (coords[0][1] - uy) ** 2
                d_last = (coords[-1][0] - ux) ** 2 + (coords[-1][1] - uy) ** 2
                if d_last < d_first:
                    coords.reverse()
            else:
                coords = [(ux, uy), (vx, vy)]
            c = np.asarray(coords, dtype=float)
            lon1.append(c[:-1, 0]); lat1.append(c[:-1, 1])
            lon2.append(c[1:, 0]); lat2.append(c[1:, 1])
            edge_idx.append(np.full(len(c) - 1, k))
            edge_tt.append(float(data["travel_time"]))

        lon1, lat1 = np.concatenate(lon1), np.concatenate(lat1)
        lon2, lat2 = np.concatenate(lon2), np.concatenate(lat2)
        edge_idx = np.concatenate(edge_idx)
        edge_tt = np.asarray(edge_tt)

        az, _, length = self.geod.inv(lon1, lat1, lon2, lat2)
        keep = length > 0.0
        lon1, lat1, lon2, lat2 = lon1[keep], lat1[keep], lon2[keep], lat2[keep]
        az, length, edge_idx = az[keep], length[keep], edge_idx[keep]

        edge_len = np.bincount(edge_idx, weights=length, minlength=len(edge_tt))
        duration = edge_tt[edge_idx] * length / edge_len[edge_idx]

        return _Polyline(lon1, lat1, lon2, lat2, length, np.mod(az, 360.0), duration)

    # -------------------------------------------------------------------------
    # Solar computation
    # -------------------------------------------------------------------------

    def _compute_segments(self, poly: _Polyline, depart: pd.Timestamp, scale: float) -> pd.DataFrame:
        """Vectorized solar geometry and exposure for every segment, evaluated at its midpoint."""
        dt_s = poly.duration_s * scale
        t_mid_s = np.cumsum(dt_s) - dt_s / 2.0
        times = pd.DatetimeIndex(depart + pd.to_timedelta(t_mid_s, unit="s"))

        lat_mid = (poly.lat1 + poly.lat2) / 2.0
        lon_mid = (poly.lon1 + poly.lon2) / 2.0

        sol = pvlib.solarposition.get_solarposition(times, lat_mid, lon_mid)
        el = sol["apparent_elevation"].to_numpy()
        az = sol["azimuth"].to_numpy()
        zen = sol["apparent_zenith"].to_numpy()

        dni = self._direct_normal_irradiance(times, zen, lat_mid, lon_mid)
        up = el > 0.0
        dni = np.where(up, dni, 0.0)

        # cos(incidence) on the left vertical window; right window is the negative.
        cos_left = np.where(up, np.cos(np.radians(el)) * np.sin(np.radians(poly.heading_deg - az)), 0.0)
        exp_left = dni * np.clip(cos_left, 0.0, None) * dt_s / 3600.0
        exp_right = dni * np.clip(-cos_left, 0.0, None) * dt_s / 3600.0

        side = np.where(
            (dni <= 0.0) | (cos_left == 0.0), "none",
            np.where(cos_left > 0.0, "right", "left"),   # sun on the left -> sit right
        )

        return pd.DataFrame({
            "time": times,
            "lat1": poly.lat1, "lon1": poly.lon1,
            "lat2": poly.lat2, "lon2": poly.lon2,
            "length_m": poly.length_m,
            "cum_m": np.cumsum(poly.length_m),
            "duration_s": dt_s,
            "heading_deg": poly.heading_deg,
            "sun_azimuth_deg": az,
            "sun_elevation_deg": el,
            "dni_w_m2": dni,
            "cos_incidence_left": cos_left,
            "exposure_left_wh_m2": exp_left,
            "exposure_right_wh_m2": exp_right,
            "side_less_sun": side,
        })

    def _direct_normal_irradiance(
        self, times: pd.DatetimeIndex, apparent_zenith: np.ndarray,
        lat: np.ndarray, lon: np.ndarray,
    ) -> np.ndarray:
        """Clear-sky DNI (Ineichen-Perez) or a constant 1.0 in geometric mode.

        Linke turbidity is looked up once at the route centroid. Its spatial resolution
        (1/12°) makes per-segment lookup unnecessary for regional routes.
        """
        if self.weighting == "geometric":
            return np.ones(len(times))

        am_rel = pvlib.atmosphere.get_relative_airmass(apparent_zenith)
        am_abs = pvlib.atmosphere.get_absolute_airmass(am_rel)
        try:
            tl = pvlib.clearsky.lookup_linke_turbidity(
                times, float(np.mean(lat)), float(np.mean(lon))
            ).to_numpy()
        except Exception as exc:  # missing h5py or data file
            logger.warning("Linke turbidity lookup failed (%s); using TL=3.0.", exc)
            tl = np.full(len(times), 3.0)
        dni_extra = pvlib.irradiance.get_extra_radiation(times).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):  # sun below horizon -> NaN, zeroed below
            cs = pvlib.clearsky.ineichen(apparent_zenith, am_abs, tl, dni_extra=dni_extra)
        return np.nan_to_num(np.asarray(cs["dni"], dtype=float), nan=0.0)

    def _summarize(
        self,
        df: pd.DataFrame,
        places: list[ResolvedPlace],
        depart: pd.Timestamp,
        arrive: pd.Timestamp,
        routed_s: float,
        scale: float,
    ) -> dict[str, object]:
        L = float(df["exposure_left_wh_m2"].sum())
        R = float(df["exposure_right_wh_m2"].sum())
        total_t = float(df["duration_s"].sum())

        # balance > 0: right window gets more sun -> sit left.
        balance = (R - L) / (R + L) if (R + L) > 0 else 0.0
        if R + L == 0:
            overall = "none"
        elif abs(balance) < self.balance_threshold:
            overall = "either"
        else:
            overall = "left" if balance > 0 else "right"

        def tshare(s: str) -> float:
            return float(df.loc[df["side_less_sun"] == s, "duration_s"].sum() / total_t)

        return {
            "origin": places[0].label,
            "destination": places[-1].label,
            "vias": [p.label for p in places[1:-1]],
            "total_distance_m": float(df["length_m"].sum()),
            "segments": int(len(df)),
            "depart_time": str(depart),
            "arrive_time": str(arrive),
            "osm_duration_min": routed_s / 60.0,
            "timing_scale": scale,
            "corridor_buffer_m": self.corridor_buffer_m,
            "weighting": self.weighting,
            "exposure_left_wh_m2": L,
            "exposure_right_wh_m2": R,
            "balance": balance,
            "time_share_sit_left": tshare("left"),
            "time_share_sit_right": tshare("right"),
            "time_share_no_direct_sun": tshare("none"),
            "overall_best_side": overall,
            "route_map_html": "",
        }


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# CLI
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate which vehicle side receives less direct sunlight along a driving route.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    place_help = "City name ('Juiz de Fora, MG, Brasil') or 'lat,lon'."
    parser.add_argument("--origin", required=True, help=place_help)
    parser.add_argument("--destination", required=True, help=place_help)
    parser.add_argument("--via", action="append", default=[], help=place_help + " Repeatable.")
    parser.add_argument("--depart", required=True, help="Departure time, e.g. '2026-02-15 14:30'.")
    parser.add_argument("--arrive", default=None,
                        help="Arrival time. If omitted, OSM travel-time estimate is used.")
    parser.add_argument("--tz", default="America/Sao_Paulo", help="Timezone for naive timestamps.")
    parser.add_argument("--corridor-buffer-m", type=float, default=8000.0,
                        help="Corridor buffer on each side of the waypoint polyline (m).")
    parser.add_argument("--weight", default="travel_time", choices=["travel_time", "length"],
                        help="Edge weight for the router.")
    parser.add_argument("--weighting", default="clearsky", choices=["clearsky", "geometric"],
                        help="Weight incidence by clear-sky DNI, or use pure geometry.")
    parser.add_argument("--balance-threshold", type=float, default=0.15,
                        help="Below this |balance| the result is 'either'.")
    parser.add_argument("--html-out", default="route.html", help="Output HTML map path.")
    parser.add_argument("--csv-out", default=None, help="Optional per-segment CSV path.")
    parser.add_argument("--no-map", action="store_true", help="Skip the HTML map.")
    parser.add_argument("--show-head", action="store_true", help="Print the first segment rows.")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="-v INFO, -vv DEBUG.")
    return parser


_PLACE_FLAGS = ("--origin", "--destination", "--via")


def _glue_negative_coords(argv: list[str]) -> list[str]:
    """Rewrite '--origin -21.7,-43.0' as '--origin=-21.7,-43.0'.

    argparse treats any token starting with '-' as an option unless it matches its
    negative-number pattern. 'lat,lon' with a comma never matches, so southern or
    western coordinates would otherwise fail with 'expected one argument'.
    """
    out: list[str] = []
    it = iter(range(len(argv)))
    for i in it:
        tok = argv[i]
        nxt = argv[i + 1] if i + 1 < len(argv) else None
        if tok in _PLACE_FLAGS and nxt is not None and _LATLON_RE.match(nxt):
            out.append(f"{tok}={nxt}")
            next(it, None)
        else:
            out.append(tok)
    return out


def _configure_logging(verbosity: int) -> None:
    level = {0: logging.WARNING, 1: logging.INFO}.get(verbosity, logging.DEBUG)
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                        datefmt="%H:%M:%S")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    args = _build_arg_parser().parse_args(_glue_negative_coords(argv))
    _configure_logging(args.verbose)

    analyzer = SunSideRouteAnalyzer(
        tz=args.tz,
        weight=args.weight,
        corridor_buffer_m=args.corridor_buffer_m,
        weighting=args.weighting,
        balance_threshold=args.balance_threshold,
    )
    try:
        df, s = analyzer.analyze(
            origin=args.origin, destination=args.destination, vias=args.via,
            depart_time=args.depart, arrive_time=args.arrive,
            save_html_map=not args.no_map, html_path=args.html_out,
        )
    except (ValueError, RuntimeError) as exc:
        logger.error("Analysis failed: %s", exc)
        return 1

    if args.csv_out:
        df.to_csv(Path(args.csv_out), index=False)

    verdict = {
        "none": "no direct sun on the side windows",
        "either": "both sides get similar sun; no clear winner",
        "left": "sit on the LEFT side",
        "right": "sit on the RIGHT side",
    }[s["overall_best_side"]]

    print(f"From:    {s['origin']}")
    print(f"To:      {s['destination']}")
    print(f"Window:  {s['depart_time']} -> {s['arrive_time']} "
          f"(OSM estimate {s['osm_duration_min']:.0f} min, scale {s['timing_scale']:.2f})")
    print(f"Route:   {s['total_distance_m'] / 1000:.1f} km, {s['segments']} segments")
    unit = "Wh/m²" if s["weighting"] == "clearsky" else "unitless"
    print(f"Exposure left/right: {s['exposure_left_wh_m2']:.0f} / {s['exposure_right_wh_m2']:.0f} {unit}"
          f"  (balance {s['balance']:+.2f})")
    print(f"Recommendation: {verdict}")
    if s["route_map_html"]:
        print(f"Map: {s['route_map_html']}")
    if args.show_head and not df.empty:
        print("\n" + df.head().to_string(index=False))
    return 0


if __name__ == "__main__":
    if os.getenv("CLEAR_SCREEN", "0") == "1":
        os.system("cls" if os.name == "nt" else "clear")
    sys.exit(main())
