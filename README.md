# SunSide Route Analyzer

SunSide Route Analyzer estimates, **along a driving route**, which side of a vehicle (**left** or **right**, relative to the direction of travel) receives **less direct sunlight** at a given date/time window.

It:
- builds a route from OpenStreetMap data (OSMnx),
- samples the route every *N meters*,
- computes the Sun position at each sampled segment (pvlib),
- compares the Sun’s horizontal direction with the vehicle’s left/right side normals,
- returns a per-segment decision plus aggregated statistics.

> This project models **direct geometric sun exposure** only. It does **not** account for shading from terrain, buildings, trees, tunnels, etc.

---

## Features

- **Segment-based analysis** (e.g., every 500 m).
- Uses **departure + arrival time** to estimate Sun position along the trip timeline.
- Supports **via points** to bias route selection:
  - Via points define a **corridor** around the intended path to avoid undesired alternative routes.
- Generates an **interactive HTML map** of the chosen route (`route.html`).
- Returns:
  - a `pandas.DataFrame` with per-segment results
  - a summary dict with route-level statistics

---

## Requirements

- Python 3.10+ (recommended: 3.11/3.12)
- Packages:
  - `osmnx`, `networkx`
  - `geopandas`, `shapely`, `pyproj`
  - `pvlib`
  - `folium`
  - `numpy`, `pandas`

### Install

```bash
pip install -r requirements.txt
```

### Quick Start

1) Run the example script

```bash
python main.py
```

By default, the script will:

- compute the route,
- generate route.html,
- print a summary and the first rows of the segment table.

Open route.html in your browser to inspect the route.

### Usage

The main entry point is:

```bash
df, summary = best_side_less_sun(
    origin_latlon=(LAT1, LON1),
    dest_latlon=(LAT2, LON2),
    via_latlon=[(LATv1, LONv1), (LATv2, LONv2)],  # optional
    depart_time="YYYY-MM-DD HH:MM",
    arrive_time="YYYY-MM-DD HH:MM",
    step_m=500.0,
    corridor_buffer_m=8000.0,
    save_html_map=True,
    html_path="route.html",
)
```

### Parameters

- origin_latlon, dest_latlon: (lat, lon)
- via_latlon (optional): list of intermediate (lat, lon) waypoints
  - Used to define a corridor around the intended path (reduces route alternatives).
- depart_time, arrive_time: timestamps (string or pd.Timestamp)
- tz: timezone (default: "America/Sao_Paulo")
- step_m: sampling step in meters (e.g., 500.0)
- corridor_buffer_m: corridor half-width in meters (e.g., 8000.0)
- save_html_map: writes an interactive route map
- html_path: output file for the map

### Output

`df (per-segment DataFrame)`

Columns include:

- time: estimated timestamp at segment start
- lat, lon: segment start coordinate
- segment_m, cum_m: segment length and cumulative distance
- road_bearing_deg: segment bearing (WGS84)
- sun_azimuth_deg, sun_elevation_deg: Sun position
- I_left, I_right: lateral incidence scores (higher = more sun on that side)
- side_less_sun: "left", "right", or "none"

`summary (aggregate stats)`

Includes:

- total distance and number of segments
- average effective incidence for left/right
- distance share where left/right had less sun
- overall recommended side (overall_best_side)
- map path (route_map_html)

### Interpretation Notes

- "left" / "right" are relative to the direction of travel.
- "none" means:
  - Sun below the horizon (elevation <= 0), or
  - Sun is approximately aligned with the direction of motion (front/back), producing minimal lateral incidence.

### Limitations

- No shading/occlusion modeling (terrain/buildings/trees).
- Route selection depends on OSM data availability and corridor settings.
- Assumes trip timeline advances linearly with distance (no live traffic profile).
