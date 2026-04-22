#!/bin/bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
"${PYTHON_BIN}" install.py "$@"
