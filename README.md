# grav-anom-tool_1.0
# Gravity Δg Time-Series (ICGEM, config-driven)

Minimal Linux-friendly, open-source Python tool to download **unfiltered 96×96 GRACE/GRACE‑FO** spherical harmonic fields from **ICGEM**, compute **gravity anomaly time‑series (Δg, µGal)** at a chosen point (e.g., an earthquake epicenter), and fit simple models (polynomial trend, Fourier seasonal terms, moving average).

- **Config-driven** via `config.ini` (location, dates, model settings).
- Works with **Python 3.9+**. Uses **numpy**; `pandas`/`matplotlib` are optional.
- **No build tools.** Just run the script.

> Data source directories are the GFZ RL06 unfiltered 96×96 sets for GRACE and GRACE-FO on ICGEM.

## Quick start

```bash
cd grav-anom-tool
cp config.example.ini config.ini
# edit config.ini for your epicenter, dates, model settings
python3 main.py --config config.ini
```

Outputs go into `outputs/latXX_lonYY/`: a CSV (or stdlib CSV) and plots if `matplotlib` is available.

## Config file (`config.ini`)

- `[general]`
  - `lat`, `lon`: decimal degrees (epicenter).
  - `start_date`, `end_date`: `YYYY-MM-DD` (optional).
  - `nmax`: maximum SH degree/order (default 96).
  - `out_dir`, `cache_dir`: folders for results and cached `.gfc` files.
- `[models]`
  - `polynomial_degree`: trend polynomial degree.
  - `fourier_terms`: # of sine/cosine pairs.
  - `moving_average_window`: odd integer window for smoothing.
- `[network]`
  - timeouts and retry knobs for ICGEM downloads.

## Dependencies

- **Required:** `numpy`
- **Optional:** `matplotlib` (plots), `pandas` (CSV writing convenience)

The tool avoids non-stdlib HTML parsers and ML libs to keep the footprint light.

## Notes

- Values are mean-removed µGal (`Δg - mean(Δg)`), suitable for model comparison.
- Plots are skipped automatically if `matplotlib` is not installed.
- CSV writing falls back to the stdlib `csv` module if `pandas` is not installed.
- The tool caches plain-text `.gfc` files and handles common ICGEM "HTML instead of file" cases.

## License

see `LICENSE`.


## New: Per-model anomaly report (last 3)
- Add optional keys under `[model]` in `config.ini`:
  - `wanted_poly_degrees = 1,2,3`
  - `wanted_fourier_terms = 1,2,4`
  - `wanted_ma_windows = 3,5`
  - `anomaly_sigma_threshold = 2.0`  (residuals > k·σ)
  - `take_last_anomalies = 3`
The pipeline will create `<lat_lon>_last_anomalies.csv` next to the time-series CSV.
