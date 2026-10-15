"""``python -m tin_engine``: the CLI and its plugins, as the MCP server runs it per call.

``sys.argv[0]`` is renamed so the run record's ``command`` field reads ``rasputin mesh ...``.
"""

import sys

from tin_engine.plugins import main

if __name__ == "__main__":
    sys.argv[0] = "rasputin"
    main()
