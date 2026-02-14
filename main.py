# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# SUNSIDE ROUTE ANALYZER
# REPOSITORY: https://github.com/cassianorcarneiro/sunside-route-analyzer
# CASSIANO RIBEIRO CARNEIRO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import folium
import geopandas as gpd
import networkx as nx
import osmnx as ox
import pvlib

from pyproj import Geod
from shapely.geometry import LineString, Point
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple, Union

# Types / Globals

LatLon = Tuple[float, float]
WGS84 = Geod(ellps="WGS84")

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

# Geometry / Solar helpers

def bearing_deg_wgs84(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    fwd_az, _, _ = WGS84.inv(lon1, lat1, lon2, lat2)
    return float((fwd_az + 360.0) % 360.0)  # 0=N, 90=E

def sun_horizontal_vector_EN(az_deg: float, el_deg: float) -> np.ndarray:
    """Horizontal projection of sun direction in local East-North plane."""
    if el_deg <= 0.0:
        return np.array([0.0, 0.0], dtype=float)
    az = np.deg2rad(az_deg)
    el = np.deg2rad(el_deg)
    c = np.cos(el)
    return np.array([c * np.sin(az), c * np.cos(az)], dtype=float)  # [E, N]

def build_route_linestring_projected(gdf_edges: gpd.GeoDataFrame) -> LineString:
    """Concatenate edge geometries (projected CRS) into a single LineString."""
    lines: List[LineString] = []
    for geom in gdf_edges.geometry:
        if geom is None:
            continue
        if geom.geom_type == "LineString":
            lines.append(geom)
        else:
            lines.extend(list(geom.geoms))  # MultiLineString

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

def sample_linestring_by_step_m(ls_proj: LineString, step_m: float) -> Tuple[List[Point], np.ndarray, float]:
    """
    Sample a projected LineString (meters) at fixed spacing.
    Returns: (points_proj, cumulative_distances_m, total_m)
    """
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

# OSM / Routing helpers

def graph_from_waypoint_corridor(waypoints_latlon: List[LatLon],
                                 buffer_m: float = 8000.0,
                                 network_type: str = "drive",) -> nx.MultiDiGraph:

    if len(waypoints_latlon) < 2:
        raise ValueError("Need at least origin and destination in waypoints_latlon.")

    line = LineString([(lon, lat) for (lat, lon) in waypoints_latlon])
    corridor = (
        gpd.GeoSeries([line], crs="EPSG:4326")
        .to_crs(3857)          # meters
        .buffer(buffer_m)      # corridor half-width
        .to_crs("EPSG:4326")
        .iloc[0]
    )
    return ox.graph_from_polygon(corridor, network_type=network_type)

def save_route_map_html(gdf_edges: gpd.GeoDataFrame,origin_latlon: LatLon,
                        dest_latlon: LatLon,
                        vias_latlon: Optional[List[LatLon]] = None,
                        out_html: str = "route.html",
                        zoom_start: int = 9) -> str:
    
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
        folium.CircleMarker(
            [vlat, vlon], radius=5, tooltip=f"Via {i+1}", fill=True
        ).add_to(m)

    folium.LayerControl().add_to(m)
    m.save(out_html)
    return out_html

# Main algorithm

def best_side_less_sun(origin_latlon: LatLon,
                       dest_latlon: LatLon,
                       depart_time: Union[str, pd.Timestamp],
                       arrive_time: Union[str, pd.Timestamp],
                       tz: str = "America/Sao_Paulo",
                       step_m: float = 500.0,
                       network_type: str = "drive",
                       via_latlon: Optional[List[LatLon]] = None,
                       corridor_buffer_m: float = 8000.0,
                       weight: str = "travel_time",
                       save_html_map: bool = True,
                       html_path: str = "route.html") -> Tuple[pd.DataFrame, Dict[str, Union[float, int, str]]]:

    # time handling

    depart_time = pd.Timestamp(depart_time)
    arrive_time = pd.Timestamp(arrive_time)

    if depart_time.tzinfo is None:
        depart_time = depart_time.tz_localize(tz)
    else:
        depart_time = depart_time.tz_convert(tz)

    if arrive_time.tzinfo is None:
        arrive_time = arrive_time.tz_localize(tz)
    else:
        arrive_time = arrive_time.tz_convert(tz)

    if arrive_time <= depart_time:
        raise ValueError("arrive_time must be later than depart_time")

    # settings / caching

    ox.settings.use_cache = True
    ox.settings.log_console = False

    # build corridor graph

    waypoints = [origin_latlon] + (via_latlon or []) + [dest_latlon]
    G = graph_from_waypoint_corridor(
        waypoints_latlon=waypoints,
        buffer_m=float(corridor_buffer_m),
        network_type=network_type,
    )

    G = ox.project_graph(G)  # project to meters for correct sampling/directions
    G = ox.add_edge_speeds(G)
    G = ox.add_edge_travel_times(G)

    # projected CRS

    nodes_gdf = ox.graph_to_gdfs(G, nodes=True, edges=False)
    gdf_nodes = nodes_gdf[0] if isinstance(nodes_gdf, tuple) else nodes_gdf
    crs_proj = gdf_nodes.crs

    # snap origin/destination

    (olat, olon) = origin_latlon
    (dlat, dlon) = dest_latlon
    p_orig = gpd.GeoSeries([Point(olon, olat)], crs="EPSG:4326").to_crs(crs_proj).iloc[0]
    p_dest = gpd.GeoSeries([Point(dlon, dlat)], crs="EPSG:4326").to_crs(crs_proj).iloc[0]

    n_orig = ox.nearest_nodes(G, X=p_orig.x, Y=p_orig.y)
    n_dest = ox.nearest_nodes(G, X=p_dest.x, Y=p_dest.y)

    # route

    route = nx.shortest_path(G, n_orig, n_dest, weight=weight)
    gdf_edges = ox.routing.route_to_gdf(G, route)

    if save_html_map:
        save_route_map_html(gdf_edges, origin_latlon, dest_latlon, via_latlon, out_html=html_path)

    # build a single route linestring and sample (meters)

    route_ls_proj = build_route_linestring_projected(gdf_edges)
    pts_proj, cum_m, total_m = sample_linestring_by_step_m(route_ls_proj, float(step_m))

    # lat/lon for solar position

    pts_ll = gpd.GeoSeries(pts_proj, crs=crs_proj).to_crs("EPSG:4326")

    total_dt_s = (arrive_time - depart_time).total_seconds()
    results: List[SegmentResult] = []

    for i in range(len(pts_proj) - 1):
        p1, p2 = pts_proj[i], pts_proj[i + 1]
        seg_m = float(cum_m[i + 1] - cum_m[i])
        if seg_m <= 0:
            continue

        frac = float(cum_m[i] / total_m) if total_m > 0 else 0.0
        t = depart_time + pd.to_timedelta(frac * total_dt_s, unit="s")

        # direction in projected CRS (treat as local EN for our purposes)

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
        svec = sun_horizontal_vector_EN(az, el)

        I_left = float(np.dot(svec, nL))
        I_right = float(np.dot(svec, nR))

        effL = max(0.0, I_left)
        effR = max(0.0, I_right)

        if el <= 0.0 or (effL == 0.0 and effR == 0.0):
            side = "none"
        else:
            side = "left" if effL < effR else "right"

        # bearing for reporting (WGS84)

        if i + 1 < len(pts_ll):
            lat2 = float(pts_ll.iloc[i + 1].y)
            lon2 = float(pts_ll.iloc[i + 1].x)
            b = bearing_deg_wgs84(lat, lon, lat2, lon2)
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

    df = pd.DataFrame([asdict(r) for r in results])

    if df.empty:
        return df, {
            "total_distance_m": float(total_m),
            "segments": 0,
            "step_m": float(step_m),
            "depart_time": str(depart_time),
            "arrive_time": str(arrive_time),
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

    summary: Dict[str, Union[float, int, str]] = {
        "total_distance_m": float(total_m),
        "segments": int(len(df)),
        "step_m": float(step_m),
        "depart_time": str(depart_time),
        "arrive_time": str(arrive_time),
        "corridor_buffer_m": float(corridor_buffer_m),
        "avg_effective_incidence_left": avg_eff_left,
        "avg_effective_incidence_right": avg_eff_right,
        "distance_share_less_sun_left": share_left,
        "distance_share_less_sun_right": share_right,
        "distance_share_no_direct_sun": share_none,
        "overall_best_side": overall,
        "route_map_html": html_path if save_html_map else "",
    }

    return df, summary

# CLI

def main() -> None:
    
    origin: LatLon = (-21.76499940229693, -43.34904075036292)
    destination: LatLon = (-21.363765592156483, -42.479130078913386)
    
    vias: List[LatLon] = [
        (-21.529919919177065, -42.64351463191004),
        (-21.729457992953115, -43.066059553474155),
    ]

    df, summary = best_side_less_sun(origin_latlon=origin,
                                     dest_latlon=destination,
                                     via_latlon=vias,
                                     depart_time="2026-02-14 14:00",
                                     arrive_time="2026-02-14 17:00",
                                     step_m=500.0,
                                     corridor_buffer_m=8000.0,
                                     save_html_map=True,
                                     html_path="route.html")

    print()
    print(summary)
    print()
    print(df[["index", "time", "road_bearing_deg", "sun_azimuth_deg", "sun_elevation_deg", "side_less_sun"]].head())
    print()
    print("Expected segments ~", int(float(summary["total_distance_m"]) // float(summary["step_m"])))
    print(f"Route map saved to: {summary['route_map_html']}")
    print()
    print(f"Overall best side: {summary['overall_best_side']}")

if __name__ == "__main__":
    os.system("cls" if os.name == "nt" else "clear")
    main()