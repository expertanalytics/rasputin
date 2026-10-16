# The Bergen Line inputs for increment 33

**Source.** Bane NOR, *Jernbane - Banenettverk* (Geonorge metadata
`c3da3591-cded-4584-a4b1-bc61b7d1f4f2`), national GML 3.2 in EPSG:25833:
`https://nedlasting.geonorge.no/geonorge/Samferdsel/Banenettverk/GML/Samferdsel_0000_Norge_25833_Banenettverk_GML.zip`,
fetched 2026-10-09 with a generic User-Agent (`curl -A "Mozilla/5.0"`).
The zip is 12 015 862 bytes; the GML inside is dated 2025-03-06, SHA-256
`85c6bf3298c30e888a5215010076205b3e00edc3c35f0ec773e179693bfc7a72`.

**Licence.** Contains data under the Norwegian licence for Open Government
data (NLOD) distributed by Bane NOR, `http://data.norge.no/nlod/no/1.0`.
The data were changed: links selected and written as GeoJSON, and domains
buffered from them (below).

**The selection** (`prepare.py`, beside this file; it repeats
`docs/increments/33-probes/probe_line.py` and `section.py`): links in
operation (`banestatus` `I`), tunnels kept.

| file (in `../rasputin_data/banenor_banenettverk/`, not committed) | what | SHA-256 |
|---|---|---|
| `bergensbanen.geojson` | Bergensbanen, 1 340 links | `75e72529…` |
| `bergen_line3.geojson` | Bergensbanen, Randsfjordbanen, Drammenbanen, 1 499 links | `fb51551d…` |
| `section_domain.geojson` | Geilo to Ål: x 127 000-150 000, 5 km buffer, simplified 20 m; 266.9 km², 53 vertices | `d3af5872…` |
| `corridor_domain.geojson` | Hokksund to Bergen: largest part of the three lines' 5 km buffer, simplified 20 m; 4 196 km², 1 158 vertices, 3 holes | `3e543be0…` |

The probe's corridor had 1 154 vertices; this one has 1 158, with the same
area to the km². Why the four extra vertices, not looked into (this run has
shapely 2.2.0 with GEOS 3.14.1; the probe's versions were not recorded).

**To regenerate:** download and unzip the file above into
`../rasputin_data/banenor_banenettverk/`, `cd` there, and run
`<repo>/.venv/bin/python <repo>/docs/benchmarks/2026-10-09/vtol-field/bergen-line/prepare.py`
(about 10 s).
