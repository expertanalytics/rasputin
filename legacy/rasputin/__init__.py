import os
from pathlib import Path

if "RASPUTIN_DATA_DIR" not in os.environ:
    print(f"WARNING: RASPUTIN_DATA_DIR not set, using current directory.")
rasputin_data_dir = Path(os.environ.get("RASPUTIN_DATA_DIR", "."))
