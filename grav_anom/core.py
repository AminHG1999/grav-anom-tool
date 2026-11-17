#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grav-anom-tool_1.1
Gravity Δg time-series from GRACE/GRACE-FO unfiltered 96x96 ICGEM products,
driven by a simple INI config. Linux-friendly, no build systems.

Standard library + numpy required. pandas/matplotlib are OPTIONAL.
If available, CSV and plots are saved; otherwise, only TXT summaries are written.

This module intentionally avoids external HTML parsers and scikit-learn.

Author: Amin Hosseinchi
ahosseinchi78@yahoo.com
"""

from __future__ import annotations
import os, re, sys, math, gzip, time, csv, io, datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional
from urllib import request, parse, error

import numpy as np

# ----------------------------- defaults ------------------------------------

ICGEM_SOURCES = [
    # GFZ RL06
    "https://icgem.gfz-potsdam.de/sp/01_GRACE/GFZ/GFZ Release 06/96x96/unfiltered/",
    "https://icgem.gfz-potsdam.de/sp/01_GRACE/CSR/CSR%20Release%2006.3%20(GFO)/96x96/unfiltered",
]

GM = 3.986004418e14
R_EARTH = 6378136.3
G0 = GM / R_EARTH**2

UNICODE_MINUS = "\u2212"

# ----------------------------- config --------------------------------------

@dataclass
class Config:
    out_dir: Path
    cache_dir: Path
    nmax: int
    lat: float
    lon: float
    start_date: Optional[dt.date]
    end_date: Optional[dt.date]
    poly_degree: int
    fourier_terms: int
    ma_window: int
    request_timeout: int = 60
    sleep_between_requests: float = 0.2
    max_retries: int = 3
    retry_backoff: float = 2.0

def parse_date(s: str | None) -> Optional[dt.date]:
    if not s: return None
    s = s.strip()
    if not s: return None
    return dt.date.fromisoformat(s)

def load_ini(path: Path) -> Config:
    import configparser
    cp = configparser.ConfigParser()
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    cp.read(path, encoding="utf-8")

    G = cp["general"] if cp.has_section("general") else {}
    M = cp["models"] if cp.has_section("models") else {}
    N = cp["network"] if cp.has_section("network") else {}

    out_dir = Path(G.get("out_dir", "outputs")).expanduser().resolve()
    cache_dir = Path(G.get("cache_dir", "icgem_cache")).expanduser().resolve()
    nmax = int(G.get("nmax", "96"))
    lat = float(G.get("lat"))
    lon = float(G.get("lon"))
    start_date = parse_date(G.get("start_date", ""))
    end_date = parse_date(G.get("end_date", ""))

    poly_degree = int(M.get("polynomial_degree", "1"))
    fourier_terms = int(M.get("fourier_terms", "4"))
    ma_window = int(M.get("moving_average_window", "3"))

    request_timeout = int(N.get("request_timeout", "60"))
    sleep_between_requests = float(N.get("sleep_between_requests", "0.2"))
    max_retries = int(N.get("max_retries", "3"))
    retry_backoff = float(N.get("retry_backoff", "2.0"))

    return Config(
        out_dir=out_dir, cache_dir=cache_dir, nmax=nmax, lat=lat, lon=lon,
        start_date=start_date, end_date=end_date,
        poly_degree=poly_degree, fourier_terms=fourier_terms, ma_window=ma_window,
        request_timeout=request_timeout, sleep_between_requests=sleep_between_requests,
        max_retries=max_retries, retry_backoff=retry_backoff
    )

# --------------------------- utilities -------------------------------------

def _ensure_dirs(cfg: Config):
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.cache_dir.mkdir(parents=True, exist_ok=True)

def _midpoint_date(d1: dt.date, d2: dt.date) -> dt.date:
    return d1 + (d2 - d1) / 2

def date_from_filename(name: str) -> Optional[dt.date]:
    base = os.path.basename(name)
    m = re.search(r'(\d{4})(\d{2})(\d{2})[-_](\d{4})(\d{2})(\d{2})', base)
    if m:
        d1 = dt.date(int(m[1]), int(m[2]), int(m[3]))
        d2 = dt.date(int(m[4]), int(m[5]), int(m[6]))
        return _midpoint_date(d1, d2)
    m = re.search(r'(\d{4})(\d{3})[-_](\d{4})(\d{3})', base)
    if m:
        d1 = dt.datetime.strptime(f"{m[1]}{m[2]}", "%Y%j").date()
        d2 = dt.datetime.strptime(f"{m[3]}{m[4]}", "%Y%j").date()
        return _midpoint_date(d1, d2)
    m = re.search(r'(?<!\d)(\d{4})[-_]?([01]\d)(?!\d)', base)
    if m:
        y, mth = int(m[1]), int(m[2])
        if 1 <= mth <= 12:
            return dt.date(y, mth, 15)
    return None

def parse_grace_dates_from_header(header_lines: List[str]):
    for line in header_lines:
        m = re.search(r'time\s*span:\s*(\d{4}-\d{2}-\d{2})\s*(?:to|-)\s*(\d{4}-\d{2}-\d{2})', line, re.I)
        if m:
            return dt.date.fromisoformat(m[1]), dt.date.fromisoformat(m[2])
    return None, None

def fully_normalized_legendre(nmax: int, x: float) -> np.ndarray:
    P = np.zeros((nmax + 1, nmax + 1))
    P[0, 0] = 1.0
    for n in range(1, nmax + 1):
        P[n, n] = math.sqrt((2*n + 1)/(2*n)) * math.sqrt(max(0.0, 1 - x*x)) * P[n-1, n-1]
    for n in range(1, nmax + 1):
        P[n, 0] = math.sqrt(2*n + 1) * x * P[n-1, 0] / math.sqrt(2*n - 1)
    for n in range(2, nmax + 1):
        for m in range(1, n):
            P[n, m] = (math.sqrt((2*n + 1)*(2*n - 1)/((n - m)*(n + m))) * x * P[n-1, m]
                       - math.sqrt((2*n + 1)*(n + m - 1)*(n - m - 1)/((2*n - 3)*(n - m)*(n + m))) * P[n-2, m])
    return P

def compute_delta_g_point(C: np.ndarray, S: np.ndarray, lat: float, lon: float, nmax: int) -> float:
    lat_rad, lon_rad = math.radians(lat), math.radians(lon)
    x = math.cos(math.pi/2 - lat_rad)
    P = fully_normalized_legendre(nmax, float(x))
    cos_mlon = np.cos(np.arange(nmax + 1, dtype=float) * lon_rad)
    sin_mlon = np.sin(np.arange(nmax + 1, dtype=float) * lon_rad)
    dg = 0.0
    for n in range(1, nmax + 1):
        s = 0.0
        for m in range(0, n + 1):
            s += P[n, m] * (C[n, m] * cos_mlon[m] + S[n, m] * sin_mlon[m])
        dg += (n + 1) * s
    return dg * G0 * 1e8  # µGal

def _f2f(s: str) -> float:
    s = s.strip().replace('D','E').replace('d','e').replace(UNICODE_MINUS,'-')
    if s.startswith('.'): s = '0' + s
    return float(s)

GFC_RE = re.compile(
    r'^\s*(?:gfc|grcof2?|shc)\s+'
    r'(?P<n>\d+)\s+(?P<m>\d+)\s+'
    r'(?P<C>[+\-\u2212]?\d*\.?\d+(?:[DEde][+\-\u2212]?\d+)?)\s+'
    r'(?P<S>[+\-\u2212]?\d*\.?\d+(?:[DEde][+\-\u2212]?\d+)?)'
)

def parse_gfc_coeffs(text: str, nmax: int):
    C = np.zeros((nmax + 1, nmax + 1))
    S = np.zeros((nmax + 1, nmax + 1))
    header, samples = [], []
    data_started = False
    parsed = 0
    for raw in text.splitlines():
        low = raw.lstrip().lower()
        if not data_started and not low.startswith(("gfc","grcof","shc")):
            header.append(raw)
            continue
        data_started = True
        m = GFC_RE.match(raw)
        if not m:
            if len(samples) < 5 and any(tag in low for tag in ("gfc","grcof","shc")):
                samples.append(raw)
            continue
        n = int(m.group("n")); m_ = int(m.group("m"))
        if n<=nmax and m_<=nmax:
            C[n,m_] = _f2f(m.group("C"))
            S[n,m_] = _f2f(m.group("S"))
            parsed += 1
    return C,S,header,samples,parsed

# ------------------------- resilient fetch ----------------------------------

# ICGEM can rate-limit rapid automated downloads.  Keep requests polite and,
# more importantly, do not turn an already percent-encoded URL into an invalid
# double-encoded URL while retrying.
_ICGEM_MIN_REQUEST_INTERVAL = 1.0
_LAST_HTTP_REQUEST_AT = 0.0

def _encode_sp_url(u: str) -> str:
    if "/sp/" not in u:
        return u
    prefix, rest = u.split("/sp/", 1)
    # Normalise each path segment exactly once.  The old implementation used
    # quote(..., safe=""), which changed existing "%20" into "%2520" and
    # could convert a valid ICGEM URL into a 404 fallback URL.
    encoded = "/".join(parse.quote(parse.unquote(p), safe="") for p in rest.split("/"))
    return prefix + "/sp/" + encoded

def _retry_after_seconds(exc: error.HTTPError) -> Optional[float]:
    """Return Retry-After in seconds when the server provides it."""
    value = exc.headers.get("Retry-After") if exc.headers else None
    if not value:
        return None
    value = value.strip()
    try:
        return max(0.0, float(value))
    except ValueError:
        try:
            from email.utils import parsedate_to_datetime
            when = parsedate_to_datetime(value)
            if when.tzinfo is None:
                when = when.replace(tzinfo=dt.timezone.utc)
            now = dt.datetime.now(dt.timezone.utc)
            return max(0.0, (when - now).total_seconds())
        except Exception:
            return None

def _pace_http_requests(min_interval: float = _ICGEM_MIN_REQUEST_INTERVAL):
    """Enforce a small gap between ICGEM requests to avoid HTTP 429 bursts."""
    global _LAST_HTTP_REQUEST_AT
    now = time.monotonic()
    wait = min_interval - (now - _LAST_HTTP_REQUEST_AT)
    if wait > 0:
        time.sleep(wait)
    _LAST_HTTP_REQUEST_AT = time.monotonic()

def _retry_get(url: str, timeout: int, retries: int, backoff: float) -> bytes:
    """Fetch a URL with 429-aware retry logic.

    HTTP 429 is treated differently from an ordinary failure: the valid URL is
    retried (respecting Retry-After when available) instead of immediately
    falling through to alternative URL spellings.  Permanent 4xx errors such
    as 404 are not retried.
    """
    last = None
    normal_attempts = 0
    rate_limit_attempts = 0
    max_rate_limit_attempts = max(int(retries), 6)

    while normal_attempts < max(1, int(retries)):
        _pace_http_requests()
        req = request.Request(
            url,
            headers={
                "User-Agent": "grav-anom-tool_1.1/1.1 (ICGEM downloader)",
                "Accept": "text/plain,text/*;q=0.9,*/*;q=0.5",
            },
        )
        try:
            with request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except error.HTTPError as e:
            last = e

            if e.code == 429:
                rate_limit_attempts += 1
                if rate_limit_attempts >= max_rate_limit_attempts:
                    raise
                server_wait = _retry_after_seconds(e)
                if server_wait is None:
                    # A longer fallback than the old 1/2/4 s sequence is
                    # deliberate: repeated rapid retries prolong rate limits.
                    server_wait = min(60.0, max(5.0, backoff ** rate_limit_attempts))
                print(f"  request rate-limited (HTTP 429); retrying same URL in {server_wait:.1f}s...")
                time.sleep(server_wait)
                continue

            # 404/403/etc. are permanent for this exact candidate; retrying the
            # same bad URL only adds load and can trigger further rate limits.
            if 400 <= e.code < 500:
                raise

            normal_attempts += 1
            if normal_attempts >= max(1, int(retries)):
                raise
            wait = min(60.0, backoff ** (normal_attempts - 1))
            print(f"  request failed ({e}); retrying in {wait:.1f}s...")
            time.sleep(wait)

        except (error.URLError, TimeoutError, OSError) as e:
            last = e
            normal_attempts += 1
            if normal_attempts >= max(1, int(retries)):
                raise
            wait = min(60.0, backoff ** (normal_attempts - 1))
            print(f"  request failed ({e}); retrying in {wait:.1f}s...")
            time.sleep(wait)

        except Exception as e:
            last = e
            normal_attempts += 1
            if normal_attempts >= max(1, int(retries)):
                raise
            wait = min(60.0, backoff ** (normal_attempts - 1))
            print(f"  request failed ({e}); retrying in {wait:.1f}s...")
            time.sleep(wait)

    if last is not None:
        raise last
    raise RuntimeError(f"request failed: {url}")

def _looks_like_html(text: str) -> bool:
    h = text.lstrip().lower()
    return h.startswith("<!doctype") or h.startswith("<html")

def list_icgem_files(base_url: str, cfg: Config) -> List[str]:
    print(f"Scanning ICGEM: {base_url}")
    try:
        html = _retry_get(base_url, cfg.request_timeout, cfg.max_retries, cfg.retry_backoff).decode("utf-8", "ignore")
    except Exception:
        fixed = _encode_sp_url(base_url)
        if fixed != base_url:
            print(f"  retry encoded URL: {fixed}")
            html = _retry_get(fixed, cfg.request_timeout, cfg.max_retries, cfg.retry_backoff).decode("utf-8", "ignore")
        else:
            raise
    # very simple href extractor
    hrefs = re.findall(r'href=["\\\']([^"\\\']+)["\\\']', html, flags=re.I)
    urls = []
    for href in hrefs:
        href_low = href.lower()
        if href_low.endswith(".gfc") or href_low.endswith(".gfc.gz"):
            # URL-encode any spaces or control chars
            safe_href = parse.quote(href, safe="/:+%._-")
            if href.startswith("http"):
                urls.append(safe_href)
            else:
                urls.append(parse.urljoin(base_url, safe_href))
    urls = sorted(set(urls))
    print(f"  found {len(urls)} files")
    return urls

def download_text_gfc(url: str, cfg: Config) -> Optional[str]:
    tried = set()
    # The directory listing already supplies the direct .gfc/.gfc.gz URL.
    # Keep only the original and a once-normalised encoded spelling.  The old
    # "?download" fallbacks return 404 on these ICGEM directory URLs and, when
    # combined with rate limiting, made existing epochs look unavailable.
    enc = _encode_sp_url(url)
    candidates = [url]
    if enc != url:
        candidates.append(enc)
    for u in candidates:
        if u in tried: continue
        tried.add(u)
        try:
            data = _retry_get(u, cfg.request_timeout, cfg.max_retries, cfg.retry_backoff)
        except Exception:
            continue
        if data[:2] == b"\\x1f\\x8b":
            try:
                text = gzip.decompress(data).decode("utf-8", "ignore")
            except Exception:
                continue
        else:
            text = data.decode("utf-8", "ignore")
        if _looks_like_html(text):
            continue
        low = text.lower()
        if "gfc" not in low and "grcof" not in low and "shc" not in low:
            continue
        return text
    return None

def download_to_cache(url: str, cfg: Config) -> Optional[Path]:
    fname = url.split("/")[-1]
    local = cfg.cache_dir / fname
    if local.exists():
        with open(local, "rb") as f:
            chunk = f.read(2048)
        try:
            txt = chunk.decode("utf-8", errors="ignore")
            if _looks_like_html(txt):
                local.unlink(missing_ok=True)
        except Exception:
            pass
    if not local.exists():
        print(f"  downloading: {fname}")
        text = download_text_gfc(url, cfg)
        if text is None:
            print("   !! could not retrieve plain-text data; skipping")
            return None
        with open(local, "w", encoding="utf-8") as f:
            f.write(text)
        time.sleep(cfg.sleep_between_requests)
    else:
        print(f"  cached: {fname}")
    return local

# ----------------------------- pipeline -------------------------------------

def collect_all_records(cfg: Config):
    records = []
    for base in ICGEM_SOURCES:
        try:
            urls = list_icgem_files(base, cfg)
        except Exception as e:
            print(f"WARNING: list failed for {base}: {e}")
            continue
        for url in urls:
            local = download_to_cache(url, cfg)
            if local is None:
                continue
            try:
                txt = local.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                print(f"   !! cannot read {local.name}")
                continue
            C,S,header,_,parsed = parse_gfc_coeffs(txt, cfg.nmax)
            if parsed == 0 or (not np.any(C) and not np.any(S)):
                print(f"  WARNING: parsed 0 lines → {local.name}")
                continue
            d1,d2 = parse_grace_dates_from_header(header)
            date = _midpoint_date(d1,d2) if (d1 and d2) else date_from_filename(local.name)
            if date is None:
                continue
            if (cfg.start_date and date < cfg.start_date) or (cfg.end_date and date > cfg.end_date):
                continue
            rsmall = float(np.sqrt(np.nanmean((C[2:, :2]**2 + S[2:, :2]**2)))) if C.size else 0.0
            print(f"   + epoch {date}  | parsed={parsed}  | coef RMS≈{rsmall:.3e}")
            records.append((date, C, S))
    records.sort(key=lambda x: x[0])
    print(f"Total epochs collected: {len(records)}")
    return records

def time_base_30day(dates: List[dt.date]) -> List[float]:
    return [(d - dates[0]).days / 30.0 for d in dates]

def _fourier_series(t: np.ndarray, y: np.ndarray, n_terms: int) -> np.ndarray:
    if len(t)==0: return np.array([])
    Tmax = max(t) if max(t)!=0 else 1.0
    N = len(t)
    A = np.ones((N, 2*n_terms+1))
    for k in range(1, n_terms+1):
        A[:,2*k-1] = np.sin(2*np.pi*k*t/Tmax)
        A[:,2*k]   = np.cos(2*np.pi*k*t/Tmax)
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    return A @ coeffs

def _moving_average(y: np.ndarray, window: int) -> np.ndarray:
    window = max(1, int(window)|1)
    kernel = np.ones(window, dtype=float)/window
    return np.convolve(y, kernel, mode='same')

def _rmse(y: np.ndarray, fit: np.ndarray) -> float:
    y = np.asarray(y); fit = np.asarray(fit)
    return float(np.sqrt(np.mean((y-fit)**2)))

def save_csv_or_txt(dates, vals, out_path: Path):
    try:
        import pandas as pd  # optional
        df = pd.DataFrame({"date": pd.to_datetime(dates), "delta_g_ug": vals})
        df.to_csv(out_path, index=False)
        print(f"saved CSV: {out_path}")
    except Exception as e:
        # fallback to plain CSV via csv module
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["date","delta_g_ug"])
            for d,v in zip(dates,vals):
                w.writerow([d.isoformat(), f"{v:.6f}"])
        print(f"saved CSV (stdlib): {out_path}")

def plot_models(dates: List[dt.date], vals: np.ndarray, tag: str, out_dir: Path, cfg: Config):
    t = np.array(time_base_30day(dates)); y = np.array(vals)
    poly_deg = cfg.poly_degree if len(y) > cfg.poly_degree else max(0, len(y)-1)

    models = {}
    if poly_deg >= 0:
        coeffs = np.polyfit(t,y,poly_deg); models["Polynomial"] = np.polyval(coeffs, t)
    models["Fourier"] = _fourier_series(t,y,cfg.fourier_terms)
    models["Moving Average"] = _moving_average(y,cfg.ma_window)

    # Try plotting if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,4))
        plt.plot(dates,y,linewidth=1,label="Δg (mean-removed)")
        plt.grid(True); plt.title(f"Δg anomaly (µGal) – {tag}")
        plt.xlabel("Time"); plt.ylabel("µGal"); plt.legend(); plt.tight_layout()
        ts_path = out_dir / f"{tag}_timeseries.png"; plt.savefig(ts_path,dpi=200); plt.close()
        print(f"[{tag}] saved: {ts_path}")
        for name,fit in models.items():
            if fit is None or len(fit)!=len(y): continue
            res = y - fit; rmse = _rmse(y,fit)
            fig,axs = plt.subplots(2,1,figsize=(12,6),sharex=True,gridspec_kw={'height_ratios':[2,1]})
            axs[0].plot(dates,y,'b-o',label="Δg (mean-removed)")
            axs[0].plot(dates,fit,linestyle='--',marker='^',label=f"{name} fit")
            axs[0].set_ylabel("µGal"); axs[0].set_title(f"{name} fit (RMSE={rmse:.3f})")
            axs[0].legend(); axs[0].grid(True)
            axs[1].bar(dates,res,width=10); axs[1].axhline(0,color='black',linewidth=0.8)
            axs[1].set_ylabel("Residuals (µGal)"); axs[1].set_xlabel("Time")
            axs[1].set_title(f"Residuals – {name}"); axs[1].grid(True)
            fig.autofmt_xdate(); plt.tight_layout()
            fig_path = out_dir / f"{tag}_{name.lower().replace(' ','_')}.png"
            plt.savefig(fig_path,dpi=200); plt.close()
            print(f"[{tag}] {name} RMSE={rmse:.3f} → {fig_path}")
    except Exception as e:
        print(f"matplotlib not available or plotting failed: {e}")

def compute_timeseries_for_point(lat: float, lon: float, records, nmax: int):
    rows = []
    for date,C,S in records:
        rows.append((date, compute_delta_g_point(C,S,lat,lon,nmax)))
    if not rows: return [], np.array([])
    rows.sort(); dates,vals = zip(*rows)
    vals = np.asarray(vals) - np.nanmean(vals)
    return list(dates), vals

def run_pipeline(cfg: Config):
    _ensure_dirs(cfg)
    print("Collecting UNFILTERED GRACE + GRACE-FO coefficients (nmax=96) from ICGEM...")
    records = collect_all_records(cfg)
    if not records:
        raise SystemExit("No records found from ICGEM (all sources tried).")

    tag = f"lat{cfg.lat:.3f}_lon{cfg.lon:.3f}"
    pt_dir = cfg.out_dir / tag; pt_dir.mkdir(parents=True, exist_ok=True)
    dates, vals = compute_timeseries_for_point(cfg.lat,cfg.lon,records,cfg.nmax)
    if not dates:
        print(f"No valid records for {tag}")
        return
    csv_path = pt_dir / f"{tag}.csv"
    save_csv_or_txt(dates, vals, csv_path)
    plot_models(dates, vals, tag, pt_dir, cfg)
    # Also write a short summary
    with open(pt_dir / f"{tag}_summary.txt", "w", encoding="utf-8") as f:
        f.write(f"epochs={len(dates)}  start={dates[0]}  end={dates[-1]}\\n")
        f.write("values are mean-removed µGal\\n")

    # >>> New: anomalies report (non-breaking; skips silently if something fails)
    try:
        from .config_ext import load_extended_options
        from .anomaly import AnomalyConfig, write_last_anomalies_csv
        from pathlib import Path
        ext = load_extended_options(Path('config.ini'))
        anocfg = AnomalyConfig(
            poly_degrees=ext.wanted_poly_degrees,
            fourier_terms=ext.wanted_fourier_terms,
            ma_windows=ext.wanted_ma_windows,
            sigma_threshold=ext.anomaly_sigma_threshold,
            take_last_n=ext.take_last_n
        )
        anopath = pt_dir / f"{tag}_last_anomalies.csv"
        # Load timeseries from CSV we just wrote (ensures same pipeline path even if pandas missing)
        import csv as _csv
        _dates, _vals = [], []
        with open(csv_path, "r", encoding="utf-8") as _f:
            r = _csv.reader(_f)
            header = next(r, None)
            # heuristic: first two columns are date,value
            for row in r:
                if not row: 
                    continue
                _dates.append(str(row[0]))
                _vals.append(float(row[1]))
        write_last_anomalies_csv(_dates, _vals, anopath, anocfg)
        print(f"Anomalies report written: {anopath}")
    except Exception as e:
        print(f"[warn] Anomalies step skipped: {e}")
    # <<< End new
    print("Done.")
