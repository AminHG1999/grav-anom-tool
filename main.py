#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI entry for the gravity anomaly tool.
Usage:
  python3 main.py --config config.ini
"""
import argparse
from pathlib import Path
from grav_anom.core import load_ini, run_pipeline

def main():
    p = argparse.ArgumentParser(description="GRACE Δg time-series from ICGEM (config-driven)")
    p.add_argument("--config", "-c", default="config.ini", help="Path to INI config file")
    args = p.parse_args()
    cfg = load_ini(Path(args.config))
    run_pipeline(cfg)

if __name__ == "__main__":
    main()
