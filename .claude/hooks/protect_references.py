"""PreToolUse hook: block edits to benchmark references and tolerances.

These numbers are owned by the maintainer (see CLAUDE.md). To allow an
intentional change, create the file `.claude/unlock-references`, make the
edit, and delete the file again.
"""

import json
import sys
from pathlib import Path

data = json.load(sys.stdin)
project = Path(data.get("cwd", "."))
path = data.get("tool_input", {}).get("file_path", "")
protected = "benchmarks/references/" in Path(path).as_posix()
unlocked = (project / ".claude" / "unlock-references").exists()

if protected and not unlocked:
    print(
        "Blocked: benchmarks/references/ holds maintainer-owned reference values and "
        "tolerances. Propose the change in the PR description instead of editing it.",
        file=sys.stderr,
    )
    sys.exit(2)
