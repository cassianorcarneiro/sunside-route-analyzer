# run_example.py — quick-start example. Open and press F5 in VS Code.

from sunside_route_analyzer import SunSideRouteAnalyzer, LatLon

analyzer = SunSideRouteAnalyzer(
    tz="America/Sao_Paulo",
    corridor_buffer_m=8000.0,
    step_m=500.0,
)

df, summary = analyzer.analyze(
    origin_latlon=LatLon(-21.363765592156483, -42.479130078913386),
    dest_latlon=LatLon(-21.76499940229693, -43.34904075036292),
    via_latlon=[
        LatLon(-21.529919919177065, -42.64351463191004),
        LatLon(-21.729457992953115, -43.066059553474155),
    ],
    depart_time="2026-02-15 14:30",
    arrive_time="2026-02-15 17:00",
    save_html_map=True,
    html_path="route.html",
)

print(summary)

print(f"\nTotal: {summary['total_distance_m']/1000:.1f} km")
print(f"Recommended side: {summary['overall_best_side']}\n")
df.to_csv("segments.csv", index=False)