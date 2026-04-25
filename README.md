# SunSide Route Analyzer

> Estimate which side of your vehicle receives **less direct sunlight** along a driving route, segment by segment.

A Python tool that combines OpenStreetMap routing with solar position calculations to recommend the shadier side of the car for any trip — useful for long highway drives where sun glare and heat make a real difference in comfort.

---

## ☀️ How it works

The analyzer breaks your route into ~500 m segments and, for each one, computes:

1. The **vehicle's heading** at that point (from the road geometry).
2. The **sun's azimuth and elevation** at that point and time (via [`pvlib`](https://pvlib-python.readthedocs.io)).
3. The **signed projection** of the sun vector onto each side of the car.
4. A **per-segment recommendation**: `left`, `right`, or `none` (sun behind / below horizon).

Then it aggregates all segments — weighted by distance — into one overall recommendation for the whole trip.

```
              ☀️
               \
                \  sun direction
                 \
       left ←─────●─────→ right
              vehicle
              ───→ heading
```

The math is fully offline (no APIs, no keys) once OpenStreetMap data is cached.

---

## 🧰 What's under the hood

| Component | Library |
|-----------|---------|
| Driving graph from OSM | [`OSMnx`](https://osmnx.readthedocs.io) |
| Shortest path | [`NetworkX`](https://networkx.org) |
| Solar position | [`pvlib`](https://pvlib-python.readthedocs.io) |
| Geometry & projections | [`Shapely`](https://shapely.readthedocs.io), [`GeoPandas`](https://geopandas.org), [`pyproj`](https://pyproj4.github.io/pyproj/) |
| Interactive map output | [`Folium`](https://python-visualization.github.io/folium/) |
| Reverse geocoding (optional) | [`geopy`](https://geopy.readthedocs.io) (Nominatim) |

The route is computed inside a **buffered corridor** around your origin → vias → destination polyline (default 8 km wide). This keeps OSMnx from downloading entire countries and prevents the shortest-path solver from finding wild detours.

---

## 📦 Installation

Requires **Python 3.10+**.

```bash
pip install osmnx networkx pvlib geopandas shapely pyproj folium geopy pandas numpy
```

Or using `requirements.txt`:

```bash
pip install -r requirements.txt
```

> **Note:** OSMnx pulls in `geopandas`, `shapely`, `pyproj`, and `networkx` automatically, but listing them explicitly avoids version surprises.

---

## 🚀 Usage

### Command-line interface

Minimal example — a daytime drive from Bicas (MG) to Juiz de Fora (MG):

```bash
python sunside_route_analyzer.py \
    --origin -21.7197,-43.0463 \
    --destination -21.7642,-43.3496 \
    --depart "2026-02-15 14:30" \
    --arrive "2026-02-15 16:00"
```

With intermediate waypoints, CSV export, and verbose logging:

```bash
python sunside_route_analyzer.py \
    --origin -21.36,-42.48 \
    --destination -21.76,-43.35 \
    --via -21.53,-42.64 \
    --via -21.73,-43.07 \
    --depart "2026-02-15 14:30" \
    --arrive "2026-02-15 17:00" \
    --csv-out segments.csv \
    --html-out my_route.html \
    -v
```

### CLI options

| Flag | Description |
|------|-------------|
| `--origin LAT,LON` | Origin coordinate (required) |
| `--destination LAT,LON` | Destination coordinate (required) |
| `--depart TIME` | Departure timestamp, ISO 8601 (required) |
| `--arrive TIME` | Arrival timestamp, ISO 8601 (required) |
| `--via LAT,LON` | Intermediate waypoint (repeat for multiple) |
| `--tz ZONE` | Timezone for naive timestamps (default: `America/Sao_Paulo`) |
| `--corridor-buffer-m N` | Corridor width in meters (default: 8000) |
| `--step-m N` | Sampling step in meters (default: 500) |
| `--weight {travel_time,length}` | Edge weight for shortest path (default: `travel_time`) |
| `--html-out PATH` | Output HTML map path (default: `route.html`) |
| `--csv-out PATH` | Optional CSV export of all segments |
| `--no-map` | Skip rendering the HTML map |
| `--show-head` | Print first rows of the segments DataFrame |
| `--lookup-cities` | Reverse-geocode origin/destination names (network) |
| `-v`, `-vv` | Increase logging verbosity |

### Programmatic usage

```python
from sunside_route_analyzer import SunSideRouteAnalyzer, LatLon

analyzer = SunSideRouteAnalyzer(
    tz="America/Sao_Paulo",
    corridor_buffer_m=8000,
    step_m=500,
)

df, summary = analyzer.analyze(
    origin_latlon=LatLon(-21.36, -42.48),
    dest_latlon=LatLon(-21.76, -43.35),
    via_latlon=[LatLon(-21.53, -42.64)],
    depart_time="2026-02-15 14:30",
    arrive_time="2026-02-15 17:00",
    save_html_map=True,
    html_path="route.html",
)

print(f"Recommended side: {summary['overall_best_side']}")
print(df.head())
```

---

## 📊 Output

### Summary dictionary

```python
{
    "total_distance_m": 95412.3,
    "segments": 191,
    "step_m": 500.0,
    "depart_time": "2026-02-15 14:30:00-03:00",
    "arrive_time": "2026-02-15 17:00:00-03:00",
    "corridor_buffer_m": 8000.0,
    "avg_effective_incidence_left": 0.12,
    "avg_effective_incidence_right": 0.43,
    "distance_share_less_sun_left": 0.78,
    "distance_share_less_sun_right": 0.18,
    "distance_share_no_direct_sun": 0.04,
    "overall_best_side": "left",
    "route_map_html": "route.html"
}
```

### Per-segment DataFrame

| Column | Description |
|--------|-------------|
| `index` | Segment number along the route |
| `time` | Estimated time at this segment |
| `lat`, `lon` | Segment start coordinates |
| `segment_m` | Segment length (m) |
| `cum_m` | Cumulative distance from origin (m) |
| `road_bearing_deg` | Heading at this segment (0=N, 90=E) |
| `sun_azimuth_deg` | Sun azimuth |
| `sun_elevation_deg` | Sun elevation |
| `I_left`, `I_right` | Signed sun-vector projection per side |
| `side_less_sun` | Recommendation: `left`, `right`, or `none` |

### Interactive HTML map

A Folium map with the computed route, color-coded markers for origin (green) and destination (red), and circle markers for any intermediate waypoints.

---

## 🧠 Methodology details

**Heading and side normals.** The vehicle's instantaneous direction is taken from the projected (metric) coordinates of consecutive sample points. The left and right normals are perpendicular to this direction in the East-North plane.

**Sun vector.** For each sample, `pvlib.solarposition.get_solarposition` returns azimuth and elevation. These are projected onto the horizontal plane:

```
sun_E = cos(elevation) · sin(azimuth)
sun_N = cos(elevation) · cos(azimuth)
```

When the sun is below the horizon, the vector is zero and the segment is reported as `side = "none"`.

**Per-side incidence.** The signed dot product between the sun vector and each side normal yields `I_left` and `I_right`. The "shadier side" is the one with the smaller positive (or zero) projection.

**Overall recommendation.** Aggregated as a distance-weighted mean of effective (clipped at zero) per-side incidence across all segments.

---

## ⚠️ Limitations

- **Direct sun only.** The model accounts for direct solar geometry, not for clouds, terrain shadowing, vegetation, tunnels, or buildings.
- **No banking or pitch.** The vehicle is treated as moving on a horizontal plane; steep grades are ignored.
- **OSM coverage.** Quality depends on OpenStreetMap's road data in the region. Sparse areas may produce odd routings — increase `corridor_buffer_m` if needed.
- **Timing is linear.** The estimated time at each segment is interpolated linearly between departure and arrival; it does not model traffic or stops.

---

## 📁 Project layout

```
sunside-route-analyzer/
├── sunside_route_analyzer.py    # Core class + CLI
├── README.md
└── requirements.txt
```

---

## 🛣️ Roadmap

- [ ] Optional cloud-cover and terrain shading via DEM
- [ ] Heat-map overlay on the HTML output (color-coded by recommended side)
- [ ] Per-segment annotated route map (not just origin/destination markers)
- [ ] Support for round trips and multi-leg journeys with custom timing per leg
- [ ] Pip-installable package

---

> *Built for road-trippers, motion-sickness sufferers, and anyone who's ever squinted into a setting sun for three hours straight.*