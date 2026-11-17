
from __future__ import annotations
from dataclasses import dataclass
from configparser import ConfigParser
from pathlib import Path

@dataclass
class ExtendedModelOptions:
    wanted_poly_degrees: list[int]
    wanted_fourier_terms: list[int]
    wanted_ma_windows: list[int]
    anomaly_sigma_threshold: float
    take_last_n: int

def parse_list_of_ints(value: str) -> list[int]:
    if not value:
        return []
    parts = [p.strip() for p in value.replace(";", ",").split(",")]
    out = []
    for p in parts:
        if p:
            try:
                out.append(int(p))
            except ValueError:
                pass
    return out

def load_extended_options(ini_path: Path) -> ExtendedModelOptions:
    cp = ConfigParser()
    cp.read(ini_path)
    sec = 'model'
    poly = parse_list_of_ints(cp.get(sec, 'wanted_poly_degrees', fallback=''))
    fou  = parse_list_of_ints(cp.get(sec, 'wanted_fourier_terms', fallback=''))
    maw  = parse_list_of_ints(cp.get(sec, 'wanted_ma_windows', fallback=''))
    sig  = cp.getfloat(sec, 'anomaly_sigma_threshold', fallback=2.0)
    last = cp.getint(sec, 'take_last_anomalies', fallback=3)
    return ExtendedModelOptions(poly, fou, maw, sig, last)
