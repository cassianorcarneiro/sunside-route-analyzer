# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# SUNSIDE ROUTE ANALYZER
# REPOSITORY: https://github.com/cassianorcarneiro/sunside-route-analyzer
# CASSIANO RIBEIRO CARNEIRO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

"""Estimate which side of a vehicle receives less direct sunlight along a driving route.

Combines OpenStreetMap routing (via OSMnx) with solar position calculations (via pvlib)
to produce a per-segment recommendation of the shadier side of the car.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import NamedTuple

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import pvlib
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim
from pyproj import Geod
from shapely.geometry import LineString, Point


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Logging
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

logger = logging.getLogger("sunside_route_analyzer")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Data types
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class LatLon(NamedTuple):
    """A geographic coordinate in WGS84 (latitude, longitude) in decimal degrees."""

    lat: float
    lon: float


# Type aliases for clarity at call sites
LatLonLike = tuple[float, float] | LatLon
TimeLike = str | pd.Timestamp


@dataclass(frozen=True)
class SegmentResult:
    """Per-segment output for a single sampled stretch of the route.

    Attributes:
        index:               Sequential index of the segment along the route.
        time:                Timestamp at which the vehicle is estimated to be at this segment.
        lat, lon:            Coordinates of the segment start in WGS84.
        segment_m:           Length of this segment in meters.
        cum_m:               Cumulative distance from the route origin in meters.
        road_bearing_deg:    Direction of travel (0=N, 90=E) in degrees.
        sun_azimuth_deg:     Solar azimuth at this point and time.
        sun_elevation_deg:   Solar elevation at this point and time.
        I_left, I_right:     Signed projection of the sun vector on each side normal.
        side_less_sun:       Which vehicle side receives less direct sunlight: "left", "right" or "none".
    """

    index: int
    time: pd.Timestamp
    lat: float
    lon: float
    segment_m: float
    cum_m: float
    road_bearing_deg: float
    sun_azimuth_deg: float
    sun_elevation_deg: float
    I_left: float
    I_right: float
    side_less_sun: str


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Core analyzer
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class SunSideRouteAnalyzer:
    """Segment-by-segment estimate of which side of a vehicle receives LESS direct sunlight
    along a driving route, using OpenStreetMap routing (OSMnx) + pvlib solar position.

    Routing uses a corridor around [origin, vias..., destination] to reduce unwanted alternatives.
    """

    # Class-level constants
    _MIN_ENDPOINT_DISTANCE_M = 1000.0          # origin and destination must be at least this far apart
    _MAX_SAMPLE_POINTS = 2_000_000             # safety cap to prevent runaway sampling
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
        step_m: float = 500.0,
        use_osmnx_cache: bool = True,
        log_console: bool = False,
        geopy_user_agent: str = "sun-side-route-analyzer",
        geopy_min_delay_s: float = 1.1,
    ) -> None:
        if corridor_buffer_m <= 0:
            raise ValueError("corridor_buffer_m must be positive.")
        if step_m <= 0:
            raise ValueError("step_m must be positive.")

        self.tz = tz
        self.network_type = network_type
        self.weight = weight
        self.corridor_buffer_m = float(corridor_buffer_m)
        self.step_m = float(step_m)

        self.wgs84 = Geod(ellps="WGS84")

        # OSMnx settings (module-level, applied once per process)
        ox.settings.use_cache = bool(use_osmnx_cache)
        ox.settings.log_console = bool(log_console)

        # Geocoding (lazy: only used when explicitly called)
        self._geolocator = Nominatim(user_agent=geopy_user_agent, timeout=10)
        self._reverse = RateLimiter(self._geolocator.reverse, min_delay_seconds=geopy_min_delay_s)

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def analyze(
        self,
        origin_latlon: LatLonLike,
        dest_latlon: LatLonLike,
        depart_time: TimeLike,
        arrive_time: TimeLike,
        via_latlon: list[LatLonLike] | None = None,
        save_html_map: bool = True,
        html_path: str | Path = "route.html",
        zoom_start: int = 9,
    ) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
        """Run the full analysis pipeline.

        Args:
            origin_latlon:   Origin coordinate (lat, lon) in WGS84.
            dest_latlon:     Destination coordinate (lat, lon) in WGS84.
            depart_time:     Departure time (string or pandas Timestamp).
            arrive_time:     Arrival time (string or pandas Timestamp).
            via_latlon:      Optional intermediate waypoints to constrain the corridor.
            save_html_map:   If True, save an interactive HTML map of the route.
            html_path:       Output path for the HTML map.
            zoom_start:      Initial zoom level of the saved map.

        Returns:
            (segments_dataframe, summary_dict)
        """
        origin = self._coerce_latlon(origin_latlon, "origin")
        destination = self._coerce_latlon(dest_latlon, "destination")
        vias = [self._coerce_latlon(v, f"via[{i}]") for i, v in enumerate(via_latlon or [])]

        self._validate_inputs(origin, destination, depart_time, arrive_time)
        depart, arrive = self._normalize_times(depart_time, arrive_time)
        html_path = str(html_path)

        logger.info("Building corridor graph (buffer=%.0f m)...", self.corridor_buffer_m)
        waypoints = [origin, *vias, destination]
        G = self._build_corridor_graph(waypoints)

        logger.info("Projecting graph and computing edge speeds/travel times...")
        G = ox.project_graph(G)
        G = ox.add_edge_speeds(G)
        G = ox.add_edge_travel_times(G)

        # Determine the projected CRS used by OSMnx
        nodes_gdf = ox.graph_to_gdfs(G, nodes=True, edges=False)
        gdf_nodes = nodes_gdf[0] if isinstance(nodes_gdf, tuple) else nodes_gdf
        crs_proj = gdf_nodes.crs

        logger.info("Snapping endpoints and computing shortest path (weight=%s)...", self.weight)
        n_orig = self._nearest_node(G, crs_proj, origin)
        n_dest = self._nearest_node(G, crs_proj, destination)

        try:
            route = nx.shortest_path(G, n_orig, n_dest, weight=self.weight)
        except nx.NetworkXNoPath as exc:
            raise RuntimeError(
                "No drivable path found between origin and destination. "
                "Try increasing corridor_buffer_m or adjusting waypoints."
            ) from exc

        gdf_edges = ox.routing.route_to_gdf(G, route)

        if save_html_map:
            logger.info("Saving route map to %s ...", html_path)
            self.save_route_map_html(
                gdf_edges=gdf_edges,
                origin_latlon=origin,
                dest_latlon=destination,
                vias_latlon=vias,
                out_html=html_path,
                zoom_start=zoom_start,
            )

        logger.info("Sampling route geometry at %.0f m steps...", self.step_m)
        route_ls_proj = self._build_route_linestring_projected(gdf_edges)
        pts_proj, cum_m, total_m = self._sample_linestring_by_step_m(route_ls_proj, self.step_m)
        pts_ll = gpd.GeoSeries(pts_proj, crs=crs_proj).to_crs("EPSG:4326")

        logger.info("Computing solar positions for %d samples...", len(pts_proj))
        df = self._compute_segment_results(
            pts_proj=pts_proj,
            pts_ll=pts_ll,
            cum_m=cum_m,
            total_m=total_m,
            depart_time=depart,
            arrive_time=arrive,
        )

        summary = self._summarize(
            df=df,
            total_m=total_m,
            depart_time=depart,
            arrive_time=arrive,
            save_html_map=save_html_map,
            html_path=html_path,
        )

        logger.info("Analysis complete: %d segments, total %.1f km.", len(df), total_m / 1000.0)
        return df, summary

    def municipality_from_geopy(
        self, coordinate: LatLonLike, language: str = "pt-BR"
    ) -> str | None:
        """Online reverse-geocode a municipality/city name.

        This makes a network request to Nominatim. Returns None on failure
        or when no municipality-like field is present in the address.
        """
        coord = self._coerce_latlon(coordinate, "coordinate")
        for zoom in self._DEFAULT_GEOCODER_ZOOM_LEVELS:
            try:
                loc = self._reverse(
                    (coord.lat, coord.lon),
                    exactly_one=True,
                    language=language,
                    zoom=zoom,
                    addressdetails=True,
                )
            except Exception as exc:  # geopy can raise many transient errors
                logger.debug("Reverse-geocoding failed at zoom=%d: %s", zoom, exc)
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
        gdf_edges: gpd.GeoDataFrame,
        origin_latlon: LatLonLike,
        dest_latlon: LatLonLike,
        vias_latlon: list[LatLonLike] | None = None,
        out_html: str | Path = "route.html",
        zoom_start: int = 9,
    ) -> str:
        """Render and save an interactive Folium map of the computed route."""
        origin = self._coerce_latlon(origin_latlon, "origin")
        destination = self._coerce_latlon(dest_latlon, "destination")
        vias = [self._coerce_latlon(v, f"via[{i}]") for i, v in enumerate(vias_latlon or [])]

        edges_ll = (
            gdf_edges.to_crs("EPSG:4326")
            if gdf_edges.crs and gdf_edges.crs.to_string() != "EPSG:4326"
            else gdf_edges
        )

        center = [(origin.lat + destination.lat) / 2.0, (origin.lon + destination.lon) / 2.0]
        m = folium.Map(location=center, zoom_start=zoom_start)

        folium.GeoJson(edges_ll.geometry.__geo_interface__, name="route").add_to(m)
        folium.Marker([origin.lat, origin.lon], tooltip="Origin",
                      icon=folium.Icon(color="green")).add_to(m)
        folium.Marker([destination.lat, destination.lon], tooltip="Destination",
                      icon=folium.Icon(color="red")).add_to(m)

        for i, via in enumerate(vias, start=1):
            folium.CircleMarker(
                [via.lat, via.lon], radius=5, tooltip=f"Via {i}", fill=True
            ).add_to(m)

        folium.LayerControl().add_to(m)
        out_html = str(out_html)
        m.save(out_html)
        return out_html

    # -------------------------------------------------------------------------
    # Validation and coercion
    # -------------------------------------------------------------------------

    @staticmethod
    def _coerce_latlon(value: LatLonLike | None, name: str) -> LatLon:
        """Convert tuple/NamedTuple-like inputs into a validated LatLon."""
        if value is None:
            raise ValueError(f"{name} must be provided.")
        if isinstance(value, LatLon):
            lat, lon = value.lat, value.lon
        else:
            try:
                lat, lon = value  # type: ignore[misc]
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{name} must be a (lat, lon) pair.") from exc

        if not (isinstance(lat, (int, float)) and isinstance(lon, (int, float))):
            raise ValueError(f"{name} must contain numeric (lat, lon).")
        lat_f, lon_f = float(lat), float(lon)
        if not (-90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0):
            raise ValueError(f"{name} is out of range: ({lat_f}, {lon_f}).")
        return LatLon(lat_f, lon_f)

    def _validate_inputs(
        self,
        origin: LatLon,
        destination: LatLon,
        depart_time: TimeLike,
        arrive_time: TimeLike,
    ) -> None:
        """Pure (no-network) validation of the user-supplied parameters."""
        if depart_time in (None, "") or arrive_time in (None, ""):
            raise ValueError("Departure and arrival times must be provided.")

        try:
            depart = pd.Timestamp(depart_time)
            arrive = pd.Timestamp(arrive_time)
        except (ValueError, TypeError) as exc:
            raise ValueError("Invalid datetime format for depart_time or arrive_time.") from exc

        if pd.isna(depart) or pd.isna(arrive):
            raise ValueError("depart_time and arrive_time must be valid timestamps.")
        if arrive <= depart:
            raise ValueError("arrive_time must be later than depart_time.")

        # Geometric sanity check (offline, fast). Replaces fragile city-name comparison.
        _, _, dist_m = self.wgs84.inv(
            origin.lon, origin.lat, destination.lon, destination.lat
        )
        if dist_m < self._MIN_ENDPOINT_DISTANCE_M:
            raise ValueError(
                f"Origin and destination are only {dist_m:.0f} m apart "
                f"(minimum {self._MIN_ENDPOINT_DISTANCE_M:.0f} m). "
                "Routing is not meaningful at this scale."
            )

    def _normalize_times(
        self, depart_time: TimeLike, arrive_time: TimeLike
    ) -> tuple[pd.Timestamp, pd.Timestamp]:
        """Localize naive timestamps and convert tz-aware ones to the analyzer timezone."""
        depart = pd.Timestamp(depart_time)
        arrive = pd.Timestamp(arrive_time)

        depart = depart.tz_localize(self.tz) if depart.tzinfo is None else depart.tz_convert(self.tz)
        arrive = arrive.tz_localize(self.tz) if arrive.tzinfo is None else arrive.tz_convert(self.tz)

        if arrive <= depart:
            raise ValueError("arrive_time must be later than depart_time.")
        return depart, arrive

    # -------------------------------------------------------------------------
    # Graph and geometry helpers
    # -------------------------------------------------------------------------

    def _build_corridor_graph(self, waypoints: list[LatLon]) -> nx.MultiDiGraph:
        if len(waypoints) < 2:
            raise ValueError("Need at least origin and destination to build the corridor graph.")

        line = LineString([(wp.lon, wp.lat) for wp in waypoints])
        corridor = (
            gpd.GeoSeries([line], crs="EPSG:4326")
            .to_crs(3857)
            .buffer(self.corridor_buffer_m)
            .to_crs("EPSG:4326")
            .iloc[0]
        )
        return ox.graph_from_polygon(corridor, network_type=self.network_type)

    def _nearest_node(self, G: nx.MultiDiGraph, crs_proj, latlon: LatLon) -> int:
        p = (
            gpd.GeoSeries([Point(latlon.lon, latlon.lat)], crs="EPSG:4326")
            .to_crs(crs_proj)
            .iloc[0]
        )
        return int(ox.nearest_nodes(G, X=p.x, Y=p.y))

    def _bearing_deg_wgs84(self, lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        fwd_az, _, _ = self.wgs84.inv(lon1, lat1, lon2, lat2)
        return float((fwd_az + 360.0) % 360.0)

    @staticmethod
    def _sun_horizontal_vector_EN(az_deg: float, el_deg: float) -> np.ndarray:
        """Project the sun direction onto the horizontal plane in (East, North) coordinates.

        Returns the zero vector when the sun is below the horizon.
        """
        if el_deg <= 0.0:
            return np.array([0.0, 0.0], dtype=float)
        az = np.deg2rad(az_deg)
        el = np.deg2rad(el_deg)
        c = np.cos(el)
        return np.array([c * np.sin(az), c * np.cos(az)], dtype=float)

    @staticmethod
    def _build_route_linestring_projected(gdf_edges: gpd.GeoDataFrame) -> LineString:
        """Concatenate edge geometries into a single continuous LineString."""
        lines: list[LineString] = []
        for geom in gdf_edges.geometry:
            if geom is None:
                continue
            if geom.geom_type == "LineString":
                lines.append(geom)
            else:
                lines.extend(list(geom.geoms))

        if not lines:
            raise RuntimeError("Route has no edge geometries.")

        coords: list[tuple[float, float]] = []
        for k, ls in enumerate(lines):
            c = list(ls.coords)
            if k > 0 and coords and coords[-1] == c[0]:
                coords.extend(c[1:])
            else:
                coords.extend(c)
        return LineString(coords)

    @classmethod
    def _sample_linestring_by_step_m(
        cls, ls_proj: LineString, step_m: float
    ) -> tuple[list[Point], np.ndarray, float]:
        if step_m <= 0:
            raise ValueError("step_m must be > 0")

        total = float(ls_proj.length)
        if not np.isfinite(total) or total <= 0:
            raise RuntimeError("Invalid route length (check geometries/CRS).")

        n = int(np.floor(total / step_m)) + 1
        if n > cls._MAX_SAMPLE_POINTS:
            raise ValueError(
                f"Sampling would create {n} points (max {cls._MAX_SAMPLE_POINTS}); "
                "check CRS or increase step_m."
            )
        dists = np.linspace(0.0, total, n, dtype=float)
        pts = [ls_proj.interpolate(d) for d in dists]
        return pts, dists, total

    # -------------------------------------------------------------------------
    # Per-segment computation
    # -------------------------------------------------------------------------

    def _compute_segment_results(
        self,
        pts_proj: list[Point],
        pts_ll: gpd.GeoSeries,
        cum_m: np.ndarray,
        total_m: float,
        depart_time: pd.Timestamp,
        arrive_time: pd.Timestamp,
    ) -> pd.DataFrame:
        """For each segment, compute solar geometry and the recommended shadier side."""
        total_dt_s = (arrive_time - depart_time).total_seconds()
        results: list[SegmentResult] = []

        for i in range(len(pts_proj) - 1):
            p1, p2 = pts_proj[i], pts_proj[i + 1]
            seg_m = float(cum_m[i + 1] - cum_m[i])
            if seg_m <= 0:
                continue

            # Linearly interpolate timestamp along the route by cumulative distance
            frac = float(cum_m[i] / total_m) if total_m > 0 else 0.0
            t = depart_time + pd.to_timedelta(frac * total_dt_s, unit="s")

            # Local direction in projected (metric) CRS
            dx = float(p2.x - p1.x)
            dy = float(p2.y - p1.y)
            norm = (dx * dx + dy * dy) ** 0.5
            if norm == 0:
                continue
            dE, dN = dx / norm, dy / norm

            # Side normals (left = 90° CCW, right = 90° CW from heading)
            nL = np.array([-dN, dE], dtype=float)
            nR = np.array([dN, -dE], dtype=float)

            lon = float(pts_ll.iloc[i].x)
            lat = float(pts_ll.iloc[i].y)

            solpos = pvlib.solarposition.get_solarposition(t, lat, lon)
            az = float(solpos["azimuth"].iloc[0])
            el = float(solpos["elevation"].iloc[0])
            svec = self._sun_horizontal_vector_EN(az, el)

            I_left = float(np.dot(svec, nL))
            I_right = float(np.dot(svec, nR))

            # Effective irradiance is non-negative (a side cannot receive negative sun)
            effL = max(0.0, I_left)
            effR = max(0.0, I_right)

            if el <= 0.0 or (effL == 0.0 and effR == 0.0):
                side = "none"
            else:
                side = "left" if effL < effR else "right"

            # Reporting bearing computed in geographic (WGS84) coordinates
            if i + 1 < len(pts_ll):
                lat2 = float(pts_ll.iloc[i + 1].y)
                lon2 = float(pts_ll.iloc[i + 1].x)
                bearing = self._bearing_deg_wgs84(lat, lon, lat2, lon2)
            else:
                bearing = float("nan")

            results.append(
                SegmentResult(
                    index=i,
                    time=t,
                    lat=lat,
                    lon=lon,
                    segment_m=seg_m,
                    cum_m=float(cum_m[i]),
                    road_bearing_deg=bearing,
                    sun_azimuth_deg=az,
                    sun_elevation_deg=el,
                    I_left=I_left,
                    I_right=I_right,
                    side_less_sun=side,
                )
            )

        return pd.DataFrame([asdict(r) for r in results])

    def _summarize(
        self,
        df: pd.DataFrame,
        total_m: float,
        depart_time: pd.Timestamp,
        arrive_time: pd.Timestamp,
        save_html_map: bool,
        html_path: str,
    ) -> dict[str, float | int | str]:
        """Compute distance-weighted aggregates over all segments."""
        base_summary: dict[str, float | int | str] = {
            "total_distance_m": float(total_m),
            "segments": int(len(df)),
            "step_m": float(self.step_m),
            "depart_time": str(depart_time),
            "arrive_time": str(arrive_time),
            "corridor_buffer_m": float(self.corridor_buffer_m),
            "route_map_html": html_path if save_html_map else "",
        }

        if df.empty:
            base_summary["overall_best_side"] = "none"
            return base_summary

        df = df.assign(
            eff_left=df["I_left"].clip(lower=0.0),
            eff_right=df["I_right"].clip(lower=0.0),
            w=df["segment_m"],
        )
        wsum = float(df["w"].sum())

        avg_eff_left = float((df["eff_left"] * df["w"]).sum() / wsum) if wsum > 0 else 0.0
        avg_eff_right = float((df["eff_right"] * df["w"]).sum() / wsum) if wsum > 0 else 0.0

        share_left = float(df.loc[df["side_less_sun"] == "left", "w"].sum() / wsum) if wsum > 0 else 0.0
        share_right = float(df.loc[df["side_less_sun"] == "right", "w"].sum() / wsum) if wsum > 0 else 0.0
        share_none = float(df.loc[df["side_less_sun"] == "none", "w"].sum() / wsum) if wsum > 0 else 0.0

        if avg_eff_left > 0.0 or avg_eff_right > 0.0:
            overall = "left" if avg_eff_left < avg_eff_right else "right"
        else:
            overall = "none"

        base_summary.update({
            "avg_effective_incidence_left": avg_eff_left,
            "avg_effective_incidence_right": avg_eff_right,
            "distance_share_less_sun_left": share_left,
            "distance_share_less_sun_right": share_right,
            "distance_share_no_direct_sun": share_none,
            "overall_best_side": overall,
        })
        return base_summary


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# CLI
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def _parse_latlon(text: str) -> LatLon:
    """Parse a 'lat,lon' string into a LatLon."""
    try:
        lat_str, lon_str = text.split(",")
        return LatLon(float(lat_str.strip()), float(lon_str.strip()))
    except (ValueError, AttributeError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid coordinate '{text}'. Expected format: 'lat,lon' (e.g. '-21.36,-42.48')."
        ) from exc


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate which vehicle side receives less sunlight along a driving route.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--origin", type=_parse_latlon, required=True,
                        help="Origin coordinate as 'lat,lon'.")
    parser.add_argument("--destination", type=_parse_latlon, required=True,
                        help="Destination coordinate as 'lat,lon'.")
    parser.add_argument("--depart", required=True,
                        help="Departure time (ISO 8601, e.g. '2026-02-15 14:30').")
    parser.add_argument("--arrive", required=True,
                        help="Arrival time (ISO 8601, e.g. '2026-02-15 17:00').")
    parser.add_argument("--via", type=_parse_latlon, action="append", default=[],
                        help="Intermediate waypoint as 'lat,lon'. Repeat for multiple vias.")
    parser.add_argument("--tz", default="America/Sao_Paulo",
                        help="Timezone for naive timestamps.")
    parser.add_argument("--corridor-buffer-m", type=float, default=8000.0,
                        help="Corridor buffer width in meters around the path.")
    parser.add_argument("--step-m", type=float, default=500.0,
                        help="Sampling step in meters along the route.")
    parser.add_argument("--weight", default="travel_time",
                        choices=["travel_time", "length"],
                        help="Edge weight used by the shortest-path algorithm.")
    parser.add_argument("--html-out", default="route.html",
                        help="Output path for the interactive HTML map.")
    parser.add_argument("--csv-out", default=None,
                        help="If set, save the per-segment DataFrame as CSV at this path.")
    parser.add_argument("--no-map", action="store_true",
                        help="Skip rendering the HTML map.")
    parser.add_argument("--show-head", action="store_true",
                        help="Print the first rows of the segments DataFrame.")
    parser.add_argument("--lookup-cities", action="store_true",
                        help="Reverse-geocode origin/destination cities (requires network).")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase logging verbosity (-v, -vv).")
    return parser


def _configure_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)

    analyzer = SunSideRouteAnalyzer(
        tz=args.tz,
        network_type="drive",
        weight=args.weight,
        corridor_buffer_m=args.corridor_buffer_m,
        step_m=args.step_m,
    )

    try:
        df, summary = analyzer.analyze(
            origin_latlon=args.origin,
            dest_latlon=args.destination,
            via_latlon=args.via,
            depart_time=args.depart,
            arrive_time=args.arrive,
            save_html_map=not args.no_map,
            html_path=args.html_out,
        )
    except (ValueError, RuntimeError) as exc:
        logger.error("Analysis failed: %s", exc)
        return 1

    # Optional CSV export
    if args.csv_out:
        out_path = Path(args.csv_out)
        df.to_csv(out_path, index=False)
        logger.info("Saved %d segments to %s", len(df), out_path)

    # Optional reverse geocoding (network)
    origin_label, dest_label = "origin", "destination"
    if args.lookup_cities:
        origin_label = analyzer.municipality_from_geopy(args.origin) or origin_label
        dest_label = analyzer.municipality_from_geopy(args.destination) or dest_label

    overall = summary.get("overall_best_side", "none")
    overall_human = "no significant sunlight" if overall == "none" else f"keep your seat on the {overall} side"

    print(f"From {origin_label} to {dest_label}, departing {args.depart}")
    print(f"Total distance: {summary['total_distance_m'] / 1000.0:.1f} km in {summary['segments']} segments")
    print(f"Recommendation: {overall_human}")
    if summary.get("route_map_html"):
        print(f"Route map: {summary['route_map_html']}")

    if args.show_head and not df.empty:
        print()
        print(df.head().to_string(index=False))

    return 0


if __name__ == "__main__":
    if os.getenv("CLEAR_SCREEN", "0") == "1":
        os.system("cls" if os.name == "nt" else "clear")
    sys.exit(main())