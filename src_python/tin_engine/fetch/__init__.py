"""`rasputin fetch`: copy what a mesh reads from a remote source into the tile
cache. Nothing on `rasputin mesh`'s path imports this package: meshing runs
offline from the cache, and a block the cache lacks is refused with the
fetch command to run. `http.py` is the only module that touches the network."""
