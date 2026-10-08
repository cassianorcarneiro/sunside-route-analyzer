# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# SUNSIDE ROUTE ANALYZER
# REPOSITORY: https://github.com/cassianorcarneiro/sunside-route-analyzer
# CASSIANO RIBEIRO CARNEIRO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import logging

from sunside_route_analyzer import SunSideRouteAnalyzer

logging.basicConfig(level=logging.INFO)

analyzer = SunSideRouteAnalyzer(
    tz="America/Sao_Paulo",
    corridor_buffer_m=8000.0,
    weighting="clearsky",
)

df, summary = analyzer.analyze(
    origin="Laranjal, MG, Brasil",                        # names and coordinates can be mixed
    destination="Juiz de Fora, MG, Brasil",
    vias=["Leopoldina, MG, Brasil", "Bicas, MG, Brasil"],
    depart_time="2026-10-09 06:40",
    arrive_time="2026-10-09 08:45",                     # optional; omit to use OSM travel times
    save_html_map=True,
    html_path="route.html",
)

for k, v in summary.items():
    print(f"{k:28s} {v}")

df.to_csv("segments.csv", index=False)
