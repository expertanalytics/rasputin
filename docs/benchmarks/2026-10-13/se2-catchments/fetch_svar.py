"""Fetch SMHI SVAR 2022 sub-catchments (delavrinningsområden) for the four main catchments HARO 28000, 38000, 40000, 42000 from the SE1_INT, SE2_INT and NO_FI district layers of SMHI's open WFS; one GeoJSON per main catchment (EPSG:3006), paged 100 units a request (the service breaks chunked responses above about 5 MB), 3 tries a page. Run in the output folder ($RASPUTIN_DATA/sweden_smhi_svar/se2/). Licence CC BY 4.0 SE, "Källa: SMHI". Increment 39, 2026-10-12."""
import json, sys, time, urllib.parse, urllib.request
UA = "Mozilla/5.0 (generic)"
BASE = "https://opendata-view.smhi.se/SMHI_vatten/Delavrinningsomraden_2022_{d}/wfs"
NAMES = {28000: "umealven", 38000: "angermanalven", 40000: "indalsalven", 42000: "ljungan"}
PAGE = 100
def page(d, cql, start):
    q = urllib.parse.urlencode({"service": "WFS", "version": "2.0.0", "request": "GetFeature",
        "typeNames": f"SMHI_vatten:Delavrinningsomraden_2022_{d}", "outputFormat": "application/json",
        "cql_filter": cql, "count": PAGE, "startIndex": start, "sortBy": "ARO_UUID"})
    req = urllib.request.Request(BASE.format(d=d) + "?" + q, headers={"User-Agent": UA})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                raw = r.read()
            return json.loads(raw)["features"], len(raw)
        except Exception as e:
            print("retry", d, start, e, file=sys.stderr); time.sleep(2 * (k + 1))
    raise SystemExit("gave up")
def get(d, cql):
    feats, nbytes, start = [], 0, 0
    while True:
        f, n = page(d, cql, start); feats += f; nbytes += n
        if len(f) < PAGE: return feats, nbytes
        start += PAGE
t0 = time.time(); total_bytes = 0
for haro, name in NAMES.items():
    feats = []
    for d in ("SE1_INT", "SE2_INT", "NO_FI"):
        f, n = get(d, f"HARO={haro}"); feats += f; total_bytes += n
    crs = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::3006"}}
    json.dump({"type": "FeatureCollection", "crs": crs, "features": feats}, open(f"{name}_{haro}_delaro_svar2022_3006.geojson", "w"))
    print(name, haro, len(feats), "units", f"{time.time()-t0:.1f}s", flush=True)
print("bytes", total_bytes, "seconds", round(time.time()-t0, 1))
open("NOTICE.txt", "w").write("Källa: SMHI. SVAR 2022 (version 2022_1_2) delavrinningsområden, CC BY 4.0 SE, fetched from https://opendata-view.smhi.se/SMHI_vatten/ on 2026-10-12. Per-HARO files assembled from the SE1_INT, SE2_INT and NO_FI districts (a modification).\n")
