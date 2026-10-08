# ☀️ SunSide Route Analyzer

> Estimate which side of your vehicle receives less direct sunlight along a driving route, segment by segment.

A Python tool that combines OpenStreetMap routing with solar position calculations to recommend the shadier side of the car for any trip — useful for long highway drives where sun glare and heat make a real difference in comfort.

<p align="center">
  <img alt="Stack" src="https://img.shields.io/badge/Stack-OSMnx%20%2B%20pvlib-blue?style=for-the-badge">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-green?style=for-the-badge">
  <img alt="Status" src="https://img.shields.io/badge/Status-Stable-success?style=for-the-badge">
</p>

---

## 📦 How it works

1. **Resolve places.** City names (`"Juiz de Fora, MG, Brasil"`) are geocoded with Nominatim. Raw `lat,lon` also works.
2. **Route.** The OSM drivable graph is downloaded inside a corridor around the waypoints. The route passes through every via, in order.
3. **Segment.** The route is split into its native OSM vertex-to-vertex segments. No resampling.
4. **Time.** Each segment gets a timestamp from OSM edge travel times. If you give an arrival time, travel times are rescaled uniformly to match it.
5. **Sun.** Solar position and clear-sky direct normal irradiance (DNI) are computed per segment with [`pvlib`](https://pvlib-python.readthedocs.io).
6. **Integrate.** Direct-beam radiant exposure (Wh/m²) on each vertical side window is summed over the trip.

```
              ☀️
               \
                \  sun direction
                 \
       left ←─────●─────→ right
              vehicle
              ───→ heading
```

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

The route is computed inside a **buffered corridor** around your origin → vias → destination polyline (default 8 km on each side, computed in the local UTM zone). This keeps OSMnx from downloading entire countries and prevents the shortest-path solver from finding wild detours.

---

## 📋 Prerequisites

- **Python 3.10+**
- A few hundred MB of disk space for OSMnx tile cache (grows with use)
- Internet for Nominatim geocoding and for OSM downloads on first run

---

## 🚀 Quick start

### Install

```bash
pip install osmnx networkx pvlib geopandas shapely pyproj folium geopy pandas numpy
```

Or via `requirements.txt`:

```bash
pip install -r requirements.txt
```

> **Note:** OSMnx pulls in `geopandas`, `shapely`, `pyproj`, and `networkx` automatically, but listing them explicitly avoids version surprises.

### Command-line usage

With city names:

```bash
python sunside_route_analyzer.py \
    --origin "Bicas, MG, Brasil" \
    --destination "Juiz de Fora, MG, Brasil" \
    --depart "2026-02-15 14:30"
```

Without `--arrive`, the OSM travel-time estimate defines the timing.

With coordinates, vias, an explicit arrival time and CSV export:

```bash
python sunside_route_analyzer.py \
    --origin -21.36,-42.48 \
    --destination "Juiz de Fora, MG, Brasil" \
    --via "Leopoldina, MG, Brasil" \
    --via -21.73,-43.07 \
    --depart "2026-02-15 14:30" \
    --arrive "2026-02-15 17:00" \
    --csv-out segments.csv \
    -v
```

### Programmatic usage

```python
from sunside_route_analyzer import SunSideRouteAnalyzer

analyzer = SunSideRouteAnalyzer(tz="America/Sao_Paulo")

df, summary = analyzer.analyze(
    origin="Bicas, MG, Brasil",
    destination="Juiz de Fora, MG, Brasil",
    depart_time="2026-02-15 14:30",
)

print(summary["overall_best_side"], summary["balance"])
```

---

## ⚙️ CLI options

| Flag | Description |
|------|-------------|
| `--origin PLACE` | City name or `lat,lon` (required) |
| `--destination PLACE` | City name or `lat,lon` (required) |
| `--via PLACE` | Intermediate waypoint; the route passes through it (repeatable) |
| `--depart TIME` | Departure timestamp (required) |
| `--arrive TIME` | Arrival timestamp (optional; default: OSM estimate) |
| `--tz ZONE` | Timezone for naive timestamps (default: `America/Sao_Paulo`) |
| `--corridor-buffer-m N` | Corridor buffer on each side, in meters (default: 8000) |
| `--weight {travel_time,length}` | Edge weight for the router (default: `travel_time`) |
| `--weighting {clearsky,geometric}` | Weight incidence by clear-sky DNI, or pure geometry (default: `clearsky`) |
| `--balance-threshold X` | Below this \|balance\| the verdict is `either` (default: 0.15) |
| `--html-out PATH` | Output HTML map path (default: `route.html`) |
| `--csv-out PATH` | Optional CSV export of all segments |
| `--no-map` | Skip rendering the HTML map |
| `--show-head` | Print first rows of the segments DataFrame |
| `-v`, `-vv` | Increase logging verbosity |

---

## 📊 Output

### Summary dictionary

| Key | Meaning |
|-----|---------|
| `origin`, `destination`, `vias` | Resolved place labels |
| `total_distance_m`, `segments` | Route length and number of native OSM segments |
| `depart_time`, `arrive_time` | Trip window used |
| `osm_duration_min` | Unscaled OSM travel-time estimate |
| `timing_scale` | Your duration ÷ OSM estimate (warning outside 0.6–1.8) |
| `exposure_left_wh_m2`, `exposure_right_wh_m2` | Direct-beam radiant exposure per side window |
| `balance` | (R − L) / (R + L). Positive: right window sunnier, sit left |
| `time_share_sit_left/right`, `time_share_no_direct_sun` | Fraction of trip time per local verdict |
| `overall_best_side` | `left`, `right`, `either` or `none` |

### Per-segment DataFrame

| Column | Description |
|--------|-------------|
| `time` | Timestamp at the segment midpoint |
| `lat1`, `lon1`, `lat2`, `lon2` | Segment endpoints (WGS84) |
| `length_m`, `cum_m`, `duration_s` | Length, cumulative distance, scaled duration |
| `heading_deg` | Geodesic forward azimuth (0 = N, 90 = E) |
| `sun_azimuth_deg`, `sun_elevation_deg` | Solar azimuth and refraction-corrected elevation |
| `dni_w_m2` | Clear-sky DNI (1.0 in geometric mode) |
| `cos_incidence_left` | Signed cosine of incidence on the left window |
| `exposure_left_wh_m2`, `exposure_right_wh_m2` | Exposure contributed by this segment |
| `side_less_sun` | Local verdict: `left`, `right` or `none` |

### Interactive HTML map

The route is colored by the local verdict: blue = sit left, orange = sit right, gray = no direct sun.

---

## 🧠 Methodology details

**Heading.** Geodesic forward azimuth of each native OSM segment on the WGS84 ellipsoid. Edge geometries are explicitly oriented from node `u` to node `v`.

**Incidence on a side window.** Each side window is a vertical plane. With heading *b*, solar azimuth *A* and apparent elevation *e*:

```
cos(i_left)  = cos(e) · sin(b − A)
cos(i_right) = −cos(i_left)
```

**Irradiance weighting.** In `clearsky` mode each segment is weighted by Ineichen–Perez clear-sky DNI. Linke turbidity comes from pvlib's climatology at the route centroid. This captures the air-mass attenuation of a low sun.

**Exposure.** `E_side = Σ DNI · max(0, cos(i_side)) · Δt / 3600` in Wh/m². Integration is over time, not distance, because heat and glare dose accumulate with time.

**Verdict.** `balance = (E_R − E_L) / (E_R + E_L)`. If |balance| is below the threshold, the verdict is `either`.

---

## ⚠️ Limitations

- **Direct sun only.** The model accounts for direct solar geometry, not for clouds, terrain shadowing, vegetation, tunnels, or buildings.
- **No banking or pitch.** The vehicle is treated as moving on a horizontal plane; steep grades are ignored.
- **OSM coverage.** Quality depends on OpenStreetMap's road data in the region. Sparse areas may produce odd routings — increase `corridor_buffer_m` if needed.
- **Free-flow timing.** Segment times come from OSM speed limits (or per-road-class defaults), uniformly rescaled to your arrival time. Traffic and stops are not modeled.
- **Vertical windows.** Real side windows lean inward (tumblehome) and the roof blocks a high sun. Results at high solar elevation overstate side exposure.
- **City centroid.** A city name resolves to the Nominatim point for that city, usually the town center.

---

## 📁 Project structure

```
sunside-route-analyzer/
├── sunside_route_analyzer.py    # Core class + CLI
├── README.md
└── requirements.txt
```

---

## 🛣️ Roadmap

- [ ] Optional cloud-cover and terrain shading via DEM
- [x] Route map color-coded by recommended side
- [x] City-name input via geocoding
- [ ] Automatic timezone from coordinates
- [ ] Support for round trips and multi-leg journeys with custom timing per leg
- [ ] Pip-installable package

---

## 📜 License

MIT — see `LICENSE` file.

## 👤 Author

**Cassiano Ribeiro Carneiro** — [@cassianorcarneiro](https://github.com/cassianorcarneiro)

---

### 🤖 AI Assistance Disclosure

The codebase architecture, organizational structure, and stylistic formatting of this repository were refactored and optimized leveraging [Claude](https://www.anthropic.com/claude) by Anthropic. All core business logic and intellectual property remain the work of the repository authors and are governed by the project's license.

---

> *Built for road-trippers, motion-sickness sufferers, and anyone who's ever squinted into a setting sun for three hours straight.*
