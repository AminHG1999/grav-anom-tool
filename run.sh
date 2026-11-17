#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 main.py --config "${1:-config.ini}"
