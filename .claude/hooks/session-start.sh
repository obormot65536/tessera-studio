#!/bin/bash
# Install the toolchain for Claude Code cloud sessions: CalculiX + the uv workspace.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

if ! command -v ccx >/dev/null 2>&1; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq
  apt-get install -y -qq calculix-ccx >/dev/null
fi

if ! command -v uv >/dev/null 2>&1; then
  python3 -m pip install --quiet uv
fi

uv sync --quiet
