# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# SUNSIDE ROUTE ANALYZER
# REPOSITORY: https://github.com/cassianorcarneiro/sunside-route-analyzer
# CASSIANO RIBEIRO CARNEIRO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

from __future__ import annotations

import os
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple, Union

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

LatLon = Tuple[float, float]

@dataclass
class SegmentResult:
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
    side_less_sun: str  # "left" | "right" | "none"

class SunSideRouteAnalyzer:
    """
    Segment-by-segment estimate of which side of a vehicle receives LESS direct sunlight
    along a driving route, using OpenStreetMap routing (OSMnx) + pvlib solar position.

    Routing uses a corridor around [origin, vias..., destination] to reduce unwanted alternatives.
    """

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
        self.tz = tz
        self.network_type = network_type
        self.weight = weight
        self.corridor_buffer_m = float(corridor_buffer_m)
        self.step_m = float(step_m)

        self.wgs84 = Geod(ellps="WGS84")

        # OSMnx settings
        ox.settings.use_cache = bool(use_osmnx_cache)
        ox.settings.log_console = bool(log_console)

        # Geocoding (optional)
        self._geolocator = Nominatim(user_agent=geopy_user_agent, timeout=10)
        self._reverse = RateLimiter(self._geolocator.reverse, min_delay_seconds=geopy_min_delay_s)

    # Public API

    def analyze(self,
                origin_latlon: LatLon,
                dest_latlon: LatLon,
                depart_time: Union[str, pd.Timestamp],
                arrive_time: Union[str, pd.Timestamp],
                via_latlon: Optional[List[LatLon]] = None,
                save_html_map: bool = True,
                html_path: str = "route.html",
                zoom_start: int = 9) -> Tuple[pd.DataFrame, Dict[str, Union[float, int, str]]]:
        
        self._validate_inputs(origin_latlon, dest_latlon, depart_time, arrive_time)

        depart_time, arrive_time = self._normalize_times(depart_time, arrive_time)

        waypoints = [origin_latlon] + (via_latlon or []) + [dest_latlon]
        G = self._build_corridor_graph(waypoints)

        # project for metric geometry
        G = ox.project_graph(G)
        G = ox.add_edge_speeds(G)
        G = ox.add_edge_travel_times(G)

        # projected CRS
        nodes_gdf = ox.graph_to_gdfs(G, nodes=True, edges=False)
        gdf_nodes = nodes_gdf[0] if isinstance(nodes_gdf, tuple) else nodes_gdf
        crs_proj = gdf_nodes.crs

        # snap endpoints
        n_orig = self._nearest_node(G, crs_proj, origin_latlon)
        n_dest = self._nearest_node(G, crs_proj, dest_latlon)

        route = nx.shortest_path(G, n_orig, n_dest, weight=self.weight)
        gdf_edges = ox.routing.route_to_gdf(G, route)

        if save_html_map:
            self.save_route_map_html(
                gdf_edges=gdf_edges,
                origin_latlon=origin_latlon,
                dest_latlon=dest_latlon,
                vias_latlon=via_latlon,
                out_html=html_path,
                zoom_start=zoom_start,
            )

        # sample geometry in meters
        route_ls_proj = self._build_route_linestring_projected(gdf_edges)
        pts_proj, cum_m, total_m = self._sample_linestring_by_step_m(route_ls_proj, self.step_m)

        # convert sampled points to lat/lon for solar computation
        pts_ll = gpd.GeoSeries(pts_proj, crs=crs_proj).to_crs("EPSG:4326")

        df = self._compute_segment_results(
            pts_proj=pts_proj,
            pts_ll=pts_ll,
            cum_m=cum_m,
            total_m=total_m,
            depart_time=depart_time,
            arrive_time=arrive_time,
        )

        summary = self._summarize(
            df=df,
            total_m=total_m,
            depart_time=depart_time,
            arrive_time=arrive_time,
            save_html_map=save_html_map,
            html_path=html_path,
        )

        return df, summary

    def municipality_from_geopy(self, coordinate: LatLon, language: str = "pt-BR") -> Optional[str]:
        """
        Online reverse-geocode municipality/city name (best-effort, may return None).
        """
        MUNICIPAL_KEYS = (
            "city", "town", "village", "municipality",
            "city_district", "district", "borough",
            "county", "state_district",
        )

        lat, lon = coordinate
        for zoom in (10, 12, 8):
            loc = self._reverse(
                (lat, lon),
                exactly_one=True,
                language=language,
                zoom=zoom,
                addressdetails=True,
            )
            if not loc:
                continue
            addr = loc.raw.get("address", {})
            for k in MUNICIPAL_KEYS:
                if addr.get(k):
                    return addr[k]
        return None

    def save_route_map_html(
        self,
        gdf_edges: gpd.GeoDataFrame,
        origin_latlon: LatLon,
        dest_latlon: LatLon,
        vias_latlon: Optional[List[LatLon]] = None,
        out_html: str = "route.html",
        zoom_start: int = 9,
    ) -> str:
        edges_ll = (
            gdf_edges.to_crs("EPSG:4326")
            if gdf_edges.crs and gdf_edges.crs.to_string() != "EPSG:4326"
            else gdf_edges
        )

        (olat, olon) = origin_latlon
        (dlat, dlon) = dest_latlon
        m = folium.Map(location=[(olat + dlat) / 2, (olon + dlon) / 2], zoom_start=zoom_start)

        folium.GeoJson(edges_ll.geometry.__geo_interface__, name="route").add_to(m)
        folium.Marker([olat, olon], tooltip="Origin").add_to(m)
        folium.Marker([dlat, dlon], tooltip="Destination").add_to(m)

        for i, (vlat, vlon) in enumerate(vias_latlon or []):
            folium.CircleMarker([vlat, vlon], radius=5, tooltip=f"Via {i+1}", fill=True).add_to(m)

        folium.LayerControl().add_to(m)
        m.save(out_html)
        return out_html

    # Internal helpers

    def _normalize_times(self,
                         depart_time: Union[str, pd.Timestamp],
                         arrive_time: Union[str, pd.Timestamp]) -> Tuple[pd.Timestamp, pd.Timestamp]:
        
        depart = pd.Timestamp(depart_time)
        arrive = pd.Timestamp(arrive_time)

        if depart.tzinfo is None:
            depart = depart.tz_localize(self.tz)
        else:
            depart = depart.tz_convert(self.tz)

        if arrive.tzinfo is None:
            arrive = arrive.tz_localize(self.tz)
        else:
            arrive = arrive.tz_convert(self.tz)

        if arrive <= depart:
            raise ValueError("arrive_time must be later than depart_time")

        return depart, arrive

    def _build_corridor_graph(self, waypoints_latlon: List[LatLon]) -> nx.MultiDiGraph:
        if len(waypoints_latlon) < 2:
            raise ValueError("Need at least origin and destination to build the corridor graph.")

        line = LineString([(lon, lat) for (lat, lon) in waypoints_latlon])
        corridor = (
            gpd.GeoSeries([line], crs="EPSG:4326")
            .to_crs(3857)
            .buffer(self.corridor_buffer_m)
            .to_crs("EPSG:4326")
            .iloc[0]
        )
        return ox.graph_from_polygon(corridor, network_type=self.network_type)

    def _nearest_node(self, G: nx.MultiDiGraph, crs_proj, latlon: LatLon) -> int:
        lat, lon = latlon
        p = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326").to_crs(crs_proj).iloc[0]
        return int(ox.nearest_nodes(G, X=p.x, Y=p.y))

    def _bearing_deg_wgs84(self, lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        fwd_az, _, _ = self.wgs84.inv(lon1, lat1, lon2, lat2)
        return float((fwd_az + 360.0) % 360.0)

    @staticmethod
    def _sun_horizontal_vector_EN(az_deg: float, el_deg: float) -> np.ndarray:
        if el_deg <= 0.0:
            return np.array([0.0, 0.0], dtype=float)
        az = np.deg2rad(az_deg)
        el = np.deg2rad(el_deg)
        c = np.cos(el)
        return np.array([c * np.sin(az), c * np.cos(az)], dtype=float)

    @staticmethod
    def _build_route_linestring_projected(gdf_edges: gpd.GeoDataFrame) -> LineString:
        lines: List[LineString] = []
        for geom in gdf_edges.geometry:
            if geom is None:
                continue
            if geom.geom_type == "LineString":
                lines.append(geom)
            else:
                lines.extend(list(geom.geoms))

        if not lines:
            raise RuntimeError("Route has no edge geometries.")

        coords: List[Tuple[float, float]] = []
        for k, ls in enumerate(lines):
            c = list(ls.coords)
            if k > 0 and coords and coords[-1] == c[0]:
                coords.extend(c[1:])
            else:
                coords.extend(c)

        return LineString(coords)

    @staticmethod
    def _sample_linestring_by_step_m(
        ls_proj: LineString, step_m: float
    ) -> Tuple[List[Point], np.ndarray, float]:
        if step_m <= 0:
            raise ValueError("step_m must be > 0")

        total = float(ls_proj.length)
        if not np.isfinite(total) or total <= 0:
            raise RuntimeError("Invalid route length (check geometries/CRS).")

        n = int(np.floor(total / step_m)) + 1
        if n > 2_000_000:
            raise ValueError(f"Sampling would create {n} points; check CRS/step_m.")
        dists = np.linspace(0.0, total, n, dtype=float)
        pts = [ls_proj.interpolate(d) for d in dists]
        return pts, dists, total

    def _compute_segment_results(self,
                                 pts_proj: List[Point],
                                 pts_ll: gpd.GeoSeries,
                                 cum_m: np.ndarray,
                                 total_m: float,
                                 depart_time: pd.Timestamp,
                                 arrive_time: pd.Timestamp) -> pd.DataFrame:
        
        total_dt_s = (arrive_time - depart_time).total_seconds()
        results: List[SegmentResult] = []

        for i in range(len(pts_proj) - 1):
            p1, p2 = pts_proj[i], pts_proj[i + 1]
            seg_m = float(cum_m[i + 1] - cum_m[i])
            if seg_m <= 0:
                continue

            frac = float(cum_m[i] / total_m) if total_m > 0 else 0.0
            t = depart_time + pd.to_timedelta(frac * total_dt_s, unit="s")

            # local direction (projected CRS, meters)
            dx = float(p2.x - p1.x)
            dy = float(p2.y - p1.y)
            norm = (dx * dx + dy * dy) ** 0.5
            if norm == 0:
                continue
            dE, dN = dx / norm, dy / norm

            # side normals
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

            effL = max(0.0, I_left)
            effR = max(0.0, I_right)

            if el <= 0.0 or (effL == 0.0 and effR == 0.0):
                side = "none"
            else:
                side = "left" if effL < effR else "right"

            # reporting bearing in WGS84
            if i + 1 < len(pts_ll):
                lat2 = float(pts_ll.iloc[i + 1].y)
                lon2 = float(pts_ll.iloc[i + 1].x)
                b = self._bearing_deg_wgs84(lat, lon, lat2, lon2)
            else:
                b = float("nan")

            results.append(
                SegmentResult(
                    index=i,
                    time=t,
                    lat=lat,
                    lon=lon,
                    segment_m=seg_m,
                    cum_m=float(cum_m[i]),
                    road_bearing_deg=b,
                    sun_azimuth_deg=az,
                    sun_elevation_deg=el,
                    I_left=I_left,
                    I_right=I_right,
                    side_less_sun=side,
                )
            )

        return pd.DataFrame([asdict(r) for r in results])

    def _summarize(self,
                   df: pd.DataFrame,
                   total_m: float,
                   depart_time: pd.Timestamp,
                   arrive_time: pd.Timestamp,
                   save_html_map: bool,
                   html_path: str) -> Dict[str, Union[float, int, str]]:
        
        if df.empty:
            return {
                "total_distance_m": float(total_m),
                "segments": 0,
                "step_m": float(self.step_m),
                "depart_time": str(depart_time),
                "arrive_time": str(arrive_time),
                "corridor_buffer_m": float(self.corridor_buffer_m),
                "overall_best_side": "none",
                "route_map_html": html_path if save_html_map else "",
            }

        df["eff_left"] = df["I_left"].clip(lower=0.0)
        df["eff_right"] = df["I_right"].clip(lower=0.0)
        df["w"] = df["segment_m"]

        wsum = float(df["w"].sum())
        avg_eff_left = float((df["eff_left"] * df["w"]).sum() / wsum) if wsum > 0 else 0.0
        avg_eff_right = float((df["eff_right"] * df["w"]).sum() / wsum) if wsum > 0 else 0.0

        share_left = float(df.loc[df["side_less_sun"] == "left", "w"].sum() / wsum) if wsum > 0 else 0.0
        share_right = float(df.loc[df["side_less_sun"] == "right", "w"].sum() / wsum) if wsum > 0 else 0.0
        share_none = float(df.loc[df["side_less_sun"] == "none", "w"].sum() / wsum) if wsum > 0 else 0.0

        overall = "none"
        if avg_eff_left > 0.0 or avg_eff_right > 0.0:
            overall = "left" if avg_eff_left < avg_eff_right else "right"

        return {
            "total_distance_m": float(total_m),
            "segments": int(len(df)),
            "step_m": float(self.step_m),
            "depart_time": str(depart_time),
            "arrive_time": str(arrive_time),
            "corridor_buffer_m": float(self.corridor_buffer_m),
            "avg_effective_incidence_left": avg_eff_left,
            "avg_effective_incidence_right": avg_eff_right,
            "distance_share_less_sun_left": share_left,
            "distance_share_less_sun_right": share_right,
            "distance_share_no_direct_sun": share_none,
            "overall_best_side": overall,
            "route_map_html": html_path if save_html_map else "",
        }
    
    def _validate_inputs(
        self,
        origin: Optional[LatLon],
        destination: Optional[LatLon],
        depart_time: Union[str, pd.Timestamp, None],
        arrive_time: Union[str, pd.Timestamp, None],
    ) -> None:
        # Presence
        if origin is None or destination is None:
            raise ValueError("Origin and destination must be provided.")

        if not isinstance(origin, tuple) or not isinstance(destination, tuple) or len(origin) != 2 or len(destination) != 2:
            raise ValueError("Origin and destination must be (lat, lon) tuples.")

        # Coordinate range validation
        def _check_coord(name: str, coord: LatLon) -> None:
            lat, lon = coord
            if not (isinstance(lat, (int, float)) and isinstance(lon, (int, float))):
                raise ValueError(f"{name} must contain numeric (lat, lon).")
            if not (-90.0 <= float(lat) <= 90.0 and -180.0 <= float(lon) <= 180.0):
                raise ValueError(f"{name} is out of range: {coord}")

        _check_coord("origin", origin)
        _check_coord("destination", destination)

        # Time validation
        if depart_time is None or arrive_time is None:
            raise ValueError("Departure and arrival times must be provided.")

        try:
            depart = pd.Timestamp(depart_time)
            arrive = pd.Timestamp(arrive_time)
        except Exception as e:
            raise ValueError("Invalid datetime format for depart_time or arrive_time.") from e

        if arrive <= depart:
            raise ValueError("Arrival time must be later than departure time.")

        # Municipality validation (best-effort; requires network)
        origin_city = self.municipality_from_geopy(origin)
        dest_city = self.municipality_from_geopy(destination)

        if origin_city is not None and dest_city is not None:
            if origin_city.strip().casefold() == dest_city.strip().casefold():
                raise ValueError(f"Origin and destination municipalities are the same: '{origin_city}'.")

# CLI

def main() -> None:

    print(f'Find geographic coordinates at: https://www.google.com/maps/\n')

    origin: LatLon = (-21.363765592156483, -42.479130078913386)
    destination: LatLon = (-21.76499940229693, -43.34904075036292)

    depart_time = "2026-02-14 14:00-03:00"
    arrive_time = "2026-02-14 17:00-03:00"

    vias: List[LatLon] = [
        (-21.529919919177065, -42.64351463191004),
        (-21.729457992953115, -43.066059553474155),
    ]
    
    analyzer = SunSideRouteAnalyzer(tz="America/Sao_Paulo",
                                    network_type="drive",
                                    weight="travel_time",
                                    corridor_buffer_m=8000.0,
                                    step_m=500.0)
    
    df, summary = analyzer.analyze(
        origin_latlon=origin,
        dest_latlon=destination,
        via_latlon=vias,
        depart_time=depart_time,
        arrive_time=arrive_time,
        save_html_map=True,
        html_path="route.html",
    )

    origin_city = analyzer.municipality_from_geopy(origin)
    dest_city = analyzer.municipality_from_geopy(destination)

    print(f"From {origin_city} to {dest_city} at {depart_time}")
    print()
    print(f"Overall best side: {summary['overall_best_side'] if summary['overall_best_side'] != 'none' else 'No significant sunlight'}")
    print()
    print(f"Route map saved to: {summary['route_map_html']}")

    # Debug (optional)
    # print(df.head())

if __name__ == "__main__":
    os.system("cls" if os.name == "nt" else "clear")
    main()