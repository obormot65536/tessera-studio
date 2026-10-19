"""PostToolUse hook: format and lint-fix a Python file right after Claude edits it."""

import json
import subprocess
import sys

data = json.load(sys.stdin)
path = data.get("tool_input", {}).get("file_path", "")
if path.endswith(".py"):
    subprocess.run(["uv", "run", "--quiet", "ruff", "format", path], check=False)
    subprocess.run(["uv", "run", "--quiet", "ruff", "check", "--fix", "--quiet", path], check=False)
