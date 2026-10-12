"""``python -m tin_engine``: the CLI, as the MCP server's child runs it (increment 37, section 6.4).

``sys.argv[0]`` is renamed so the run record's ``command`` field reads ``rasputin mesh ...``.
"""

import sys

from tin_engine.cli import app

if __name__ == "__main__":
    sys.argv[0] = "rasputin"
    app()
