# SunSide Route Analyzer

Estimate which side of a vehicle gets less direct sunlight along a driving route, segment by segment.

![Python](https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white)
![OSMnx](https://img.shields.io/badge/OSMnx-routing-7EBC6F?logo=openstreetmap&logoColor=white)
![pvlib](https://img.shields.io/badge/pvlib-solar%20position-orange)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-stable-success)

## Overview

The tool routes a trip over the OpenStreetMap road network. It timestamps every road segment and computes the sun's position and clear-sky irradiance at that moment. It then integrates the direct-beam energy reaching each side window over the whole trip. The result is a recommendation (left, right or either), a per-segment table and an interactive map.

## How it works

```mermaid
flowchart LR
    P["Resolve places<br/>Nominatim or lat,lon"] --> R["Route in a corridor<br/>OSMnx + NetworkX"]
    R --> S["Native OSM segments"]
    S --> T["Timestamp segments<br/>OSM travel times"]
    T --> U["Sun position + DNI<br/>pvlib"]
    U --> E["Integrate exposure<br/>per side window"]
    E --> O["Summary, CSV, HTML map"]
```

1. **Resolve places.** City names are geocoded with Nominatim. Raw `lat,lon` also works.
2. **Route.** The drivable OSM graph is downloaded inside a buffered corridor around origin → vias → destination. The route visits every via, in order.
3. **Segment.** The route is split into its native vertex-to-vertex OSM segments. There is no resampling.
4. **Time.** Each segment gets a timestamp from OSM edge travel times. If an arrival time is given, all travel times are scaled uniformly to match it.
5. **Sun.** Solar position and clear-sky direct normal irradiance (DNI) are computed per segment.
6. **Integrate.** Direct-beam radiant exposure (Wh/m²) on each vertical side window is summed over the trip.

The corridor (8 km on each side by default, computed in the local UTM zone) bounds the OSM download. It also keeps the shortest-path solver from wandering into long detours.

## Tech stack

| Component | Library |
|---|---|
| OSM driving graph | [`OSMnx`](https://osmnx.readthedocs.io) |
| Shortest path | [`NetworkX`](https://networkx.org) |
| Solar position and clear-sky model | [`pvlib`](https://pvlib-python.readthedocs.io) |
| Geometry and projections | [`Shapely`](https://shapely.readthedocs.io), [`GeoPandas`](https://geopandas.org), [`pyproj`](https://pyproj4.github.io/pyproj/) |
| Interactive map | [`Folium`](https://python-visualization.github.io/folium/) |
| Geocoding | [`geopy`](https://geopy.readthedocs.io) (Nominatim) |

## Requirements

- Python 3.10+
- A few hundred MB of disk for the OSMnx cache, which grows with use
- Internet access for geocoding and for the first OSM download of each area

## Getting started

### Install

```bash
git clone https://github.com/cassianorcarneiro/sunside-route-analyzer.git
cd sunside-route-analyzer
pip install -r requirements.txt
```

OSMnx already depends on GeoPandas, Shapely, pyproj and NetworkX. `requirements.txt` lists them explicitly to pin compatible versions.

### Command line

With city names:

```bash
python sunside_route_analyzer.py \
    --origin "Bicas, MG, Brasil" \
    --destination "Juiz de Fora, MG, Brasil" \
    --depart "2026-02-15 14:30"
```

Without `--arrive`, the OSM travel-time estimate sets the timing.

With coordinates, vias, an arrival time and CSV export:

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

### Python API

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

## CLI reference

| Flag | Description |
|---|---|
| `--origin PLACE` | City name or `lat,lon` (required) |
| `--destination PLACE` | City name or `lat,lon` (required) |
| `--via PLACE` | Intermediate waypoint, visited in order (repeatable) |
| `--depart TIME` | Departure timestamp (required) |
| `--arrive TIME` | Arrival timestamp (default: OSM estimate) |
| `--tz ZONE` | Time zone for naive timestamps (default: `America/Sao_Paulo`) |
| `--corridor-buffer-m N` | Corridor half-width in meters (default: `8000`) |
| `--weight {travel_time,length}` | Edge weight for routing (default: `travel_time`) |
| `--weighting {clearsky,geometric}` | Weight incidence by clear-sky DNI, or by geometry only (default: `clearsky`) |
| `--balance-threshold X` | Below this \|balance\|, the verdict is `either` (default: `0.15`) |
| `--html-out PATH` | HTML map path (default: `route.html`) |
| `--csv-out PATH` | Optional CSV export of all segments |
| `--no-map` | Skip the HTML map |
| `--show-head` | Print the first rows of the segment table |
| `-v`, `-vv` | Increase log verbosity |

## Output

### Summary

| Key | Meaning |
|---|---|
| `origin`, `destination`, `vias` | Resolved place labels |
| `total_distance_m`, `segments` | Route length and number of OSM segments |
| `depart_time`, `arrive_time` | Trip window used |
| `osm_duration_min` | Unscaled OSM travel-time estimate |
| `timing_scale` | Your duration ÷ OSM estimate; a warning is logged outside 0.6–1.8 |
| `exposure_left_wh_m2`, `exposure_right_wh_m2` | Direct-beam exposure per side window |
| `balance` | $(E_R - E_L)/(E_R + E_L)$; positive means the right side is sunnier, so sit left |
| `time_share_sit_left`, `time_share_sit_right`, `time_share_no_direct_sun` | Share of trip time per local verdict |
| `overall_best_side` | `left`, `right`, `either` or `none` |

### Per-segment table

| Column | Description |
|---|---|
| `time` | Timestamp at the segment midpoint |
| `lat1`, `lon1`, `lat2`, `lon2` | Segment endpoints (WGS84) |
| `length_m`, `cum_m`, `duration_s` | Length, cumulative distance, scaled duration |
| `heading_deg` | Geodesic forward azimuth (0 = N, 90 = E) |
| `sun_azimuth_deg`, `sun_elevation_deg` | Solar azimuth and refraction-corrected elevation |
| `dni_w_m2` | Clear-sky DNI; `1.0` in geometric mode |
| `cos_incidence_left` | Signed cosine of incidence on the left window |
| `exposure_left_wh_m2`, `exposure_right_wh_m2` | Exposure contributed by this segment |
| `side_less_sun` | Local verdict: `left`, `right` or `none` |

### Map

The route is colored by local verdict: blue means sit left, orange means sit right, gray means no direct sun.

## Methodology

**Heading.** Each segment's heading $b$ is the geodesic forward azimuth on the WGS84 ellipsoid. Edge geometries are explicitly oriented from node `u` to node `v`.

**Incidence.** Each side window is modeled as a vertical plane. The left window's outward normal is horizontal, at azimuth $b - 90°$. For a sun at azimuth $A$ and apparent elevation $e$, the cosine of the incidence angle is

$$
\cos i_L = \cos e \,\cos\big(A - (b - 90°)\big) = \cos e \,\sin(b - A), \qquad \cos i_R = -\cos i_L .
$$

**Irradiance.** In `clearsky` mode, each segment is weighted by Ineichen–Perez clear-sky DNI. Linke turbidity comes from pvlib's climatology at the route centroid. This captures the air-mass attenuation of a low sun. In `geometric` mode, DNI is fixed at 1.

**Exposure.** Exposure is integrated over time, not distance, because heat and glare dose accumulate with time:

$$
E_s = \sum_k \mathrm{DNI}_k \,\max(0, \cos i_{s,k}) \,\frac{\Delta t_k}{3600} \quad [\mathrm{Wh/m^2}].
$$

**Verdict.** $\text{balance} = (E_R - E_L)/(E_R + E_L)$. If $|\text{balance}|$ is below the threshold, the verdict is `either`. If both exposures are zero, the verdict is `none`.

## Limitations

- **Direct sun only.** Clouds, terrain shadows, vegetation, tunnels and buildings are not modeled.
- **Flat-road assumption.** Vehicle pitch and roll are ignored.
- **Vertical windows.** Real side windows lean inward and the roof shades a high sun. Results overstate side exposure at high solar elevations.
- **Free-flow timing.** Segment times come from OSM speed limits or road-class defaults, scaled uniformly. Traffic and stops are not modeled.
- **OSM coverage.** Sparse regions can produce odd routes. Increase `--corridor-buffer-m` if needed.
- **City centroids.** A city name resolves to Nominatim's point for that city, usually the town center.
- **Nominatim usage policy.** The public instance allows at most 1 request per second. For heavy use, run your own instance or pass coordinates.

## Project structure

```
sunside-route-analyzer/
├── sunside_route_analyzer.py   # Analyzer class and CLI
├── requirements.txt
├── LICENSE
└── README.md
```

## Roadmap

- [x] Route map colored by recommended side
- [x] City-name input through geocoding
- [ ] Cloud cover and terrain shading from a DEM
- [ ] Time zone inferred from coordinates
- [ ] Round trips and multi-leg journeys with timing per leg
- [ ] Installable package on PyPI

## License

Released under the [MIT License](LICENSE).

## Author

**Cassiano Ribeiro Carneiro** · [@cassianorcarneiro](https://github.com/cassianorcarneiro)

## AI assistance

Claude (Anthropic) was used to refactor the structure and formatting of this repository. The core logic is the author's own work.
