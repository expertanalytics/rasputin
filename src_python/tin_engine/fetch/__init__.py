"""`rasputin fetch`: copy what a mesh reads from a remote source into the tile
cache (increment 23a-2, `docs/increments/23-basin-scale.md`, "The fetch step
and the tile cache"). Nothing on `rasputin mesh`'s path imports this package:
meshing is offline by rule (K7). `http.py` is the only module that touches
the network."""
