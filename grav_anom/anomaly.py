
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Anomaly detection utilities for grav-anom-tool.
- Works with numpy only.
- Detects anomalies on residuals for different model families and "orders".
- Writes a per-point CSV with the last N anomalies per model.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Dict, Iterable
import math, csv
from pathlib import Path
import numpy as np

@dataclass
class AnomalyConfig:
    # Which models + orders to run. Falls back to single values if lists are empty.
    poly_degrees: List[int]            # e.g., [1,2,3]
    fourier_terms: List[int]           # e.g., [1,2,4]
    ma_windows: List[int]              # e.g., [3,5,7]
    sigma_threshold: float = 2.0       # residuals exceeding k*sigma are anomalies
    take_last_n: int = 3               # how many latest anomalies to keep

def _to_float_dates(dates: List[str]) -> np.ndarray:
    # Dates are expected as "YYYY-MM-DD" or ISO; convert to ordinal float relative to first date
    # for scale stability.
    from datetime import datetime
    ts = np.array([datetime.fromisoformat(d).toordinal() for d in dates], dtype=float)
    t0 = ts[0]
    ts = ts - t0
    # normalize to ~years to avoid numeric issues (365.25 ~ 1 year)
    return ts / 365.25

def _design_poly(t: np.ndarray, deg: int) -> np.ndarray:
    # Columns: [1, t, t^2, ... , t^deg]
    cols = [np.ones_like(t)]
    for k in range(1, deg+1):
        cols.append(t**k)
    return np.column_stack(cols)

def _design_fourier(t: np.ndarray, m: int) -> np.ndarray:
    # base with intercept
    cols = [np.ones_like(t)]
    # one cycle per year base frequency; add harmonics up to m
    for k in range(1, m+1):
        cols.append(np.sin(2*np.pi*k*t))
        cols.append(np.cos(2*np.pi*k*t))
    return np.column_stack(cols)

def _moving_average(x: np.ndarray, w: int) -> np.ndarray:
    if w < 1 or w % 2 == 0:
        raise ValueError("moving average window must be odd, >=1")
    pad = w // 2
    # simple centered MA with edge replication
    xpad = np.pad(x, (pad, pad), mode='edge')
    kernel = np.ones(w)/w
    y = np.convolve(xpad, kernel, mode='valid')
    return y

def _fit_and_residuals(X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    # returns (yhat, resid, sigma_hat)
    # Solve least squares; add small ridge for safety if ill-conditioned.
    XtX = X.T @ X
    ridge = 1e-8 * np.eye(XtX.shape[0])
    beta = np.linalg.solve(XtX + ridge, X.T @ y)
    yhat = X @ beta
    resid = y - yhat
    dof = max(1, len(y) - X.shape[1])
    sigma_hat = float(np.sqrt(np.sum(resid**2) / dof))
    return yhat, resid, sigma_hat

def _pick_last_n_by_date(dates: np.ndarray, values: np.ndarray, mask: np.ndarray, n: int) -> List[int]:
    # indices of last n True in mask by date ordering (dates are already chronological)
    idx = [i for i in range(len(dates)) if mask[i]]
    return idx[-n:]  # already chronological, so take the last n

def detect_anomalies_per_models(
    dates: List[str], values: List[float], cfg: AnomalyConfig
) -> Dict[str, List[Dict[str, object]]]:
    """
    Returns a dict keyed by model label, each value is a list of up to cfg.take_last_n
    dicts: { 'date': str, 'value': float, 'ci_low': float, 'ci_high': float }
    where value is the residual (µGal) and CI is +/- 1.96*sigma_hat around 0.
    """
    t = _to_float_dates(dates)
    y = np.asarray(values, dtype=float)
    results: Dict[str, List[Dict[str, object]]] = {}

    # helper to format CI and filter anomalies
    def build_records(resid: np.ndarray, sigma_hat: float, label: str):
        thr = cfg.sigma_threshold * sigma_hat
        mask = np.abs(resid) > thr
        last_idx = _pick_last_n_by_date(np.arange(len(dates)), resid, mask, cfg.take_last_n)
        out = []
        ci = 1.96 * sigma_hat
        for i in last_idx:
            out.append({
                "date": dates[i],
                "value": float(resid[i]),
                "ci_low": float(-ci),
                "ci_high": float(ci)
            })
        results[label] = out

    # Polynomial family
    if not cfg.poly_degrees:
        cfg.poly_degrees = [1]
    for d in cfg.poly_degrees:
        X = _design_poly(t, d)
        _, resid, sigma = _fit_and_residuals(X, y)
        build_records(resid, sigma, f"poly_deg{d}")

    # Fourier family
    if not cfg.fourier_terms:
        cfg.fourier_terms = [1]
    for m in cfg.fourier_terms:
        X = _design_fourier(t, m)
        _, resid, sigma = _fit_and_residuals(X, y)
        build_records(resid, sigma, f"fourier_m{m}")

    # Moving average family (no regression; residual = y - MA(y))
    if not cfg.ma_windows:
        cfg.ma_windows = []
    for w in cfg.ma_windows:
        yhat = _moving_average(y, w)
        resid = y - yhat
        dof = max(1, len(y) - 1)
        sigma = float(np.sqrt(np.sum(resid**2)/dof))
        build_records(resid, sigma, f"moving_avg_w{w}")

    return results

def write_last_anomalies_csv(
    dates: List[str], values: List[float], out_csv: Path, cfg: AnomalyConfig
) -> Dict[str, List[Dict[str, object]]]:
    res = detect_anomalies_per_models(dates, values, cfg)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model", "date", "residual_value_µGal", "ci_low", "ci_high"])
        for model, rows in res.items():
            for r in rows:
                w.writerow([model, r["date"], f'{r["value"]:.6f}', f'{r["ci_low"]:.6f}', f'{r["ci_high"]:.6f}'])
    return res
