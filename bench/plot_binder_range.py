#!/usr/bin/env python3
"""
Read per-L or consolidated timeseries CSVs and plot Binder U vs ratio
for L = [48,32,24,16] on 15 ratios from 10/9 down to 0.9.
Saves: bench/summary_binder_range.csv and bench/binder_range_Ls.png
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

csv.field_size_limit(sys.maxsize)

bench = Path('bench')
bench.mkdir(exist_ok=True)

Ls = [48, 32, 24, 16]
ratios = np.linspace(10.0/9.0, 0.9, 15)
ratios = np.round(ratios, 6)

# Try per-L file first, otherwise fallback to consolidated ts_high_res.csv
per_L_paths = {L: bench / f'ts_L{L}_high_res.csv' for L in Ls}
consolidated = bench / 'ts_high_res.csv'

# Load data into dict[(L,ratio)] -> list of U
data = defaultdict(list)

# Helper to read a file
def read_file(path):
    if not path.exists():
        return 0
    rows = 0
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                L = int(r['L'])
                ratio = float(r['ratio'])
                U = float(r['U'])
            except Exception:
                continue
            key = (L, round(ratio,6))
            data[key].append(U)
            rows += 1
    return rows

# Read per-L files if present
total_rows = 0
for L,p in per_L_paths.items():
    n = read_file(p)
    total_rows += n

# If nothing read, fallback to consolidated
if total_rows == 0 and consolidated.exists():
    total_rows = read_file(consolidated)

print(f'Read {total_rows} rows total')

# Build summary list and plot
summary_rows = []
by_L = {}
for L in Ls:
    xs = []
    ys = []
    ses = []
    for r in ratios:
        key = (L, round(r,6))
        lst = data.get(key, [])
        if len(lst) == 0:
            # Try tolerance match within 1e-4
            matches = []
            for (LL, rr), vals in data.items():
                if LL==L and abs(rr - r) < 1e-4:
                    matches.extend(vals)
            lst = matches
        if len(lst) == 0:
            xs.append(r)
            ys.append(np.nan)
            ses.append(np.nan)
            summary_rows.append((L, r, 0, np.nan, np.nan))
        else:
            arr = np.array(lst)
            meanU = float(np.mean(arr))
            se = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr)>1 else 0.0
            xs.append(r)
            ys.append(meanU)
            ses.append(se)
            summary_rows.append((L, r, len(arr), meanU, se))
    by_L[L] = (np.array(xs), np.array(ys), np.array(ses))

# Write summary CSV
out_summary = bench / 'summary_binder_range.csv'
with open(out_summary, 'w') as f:
    f.write('L,ratio,n,mean_U,se_U\n')
    for row in summary_rows:
        L,r,n,meanU,se = row
        f.write(f"{L},{r:.6f},{n},{'' if np.isnan(meanU) else f'{meanU:.6f}'},{'' if np.isnan(se) else f'{se:.6f}'}\n")

print(f'Wrote {out_summary}')

# Plot
fig, ax = plt.subplots(figsize=(10,6))
colors = plt.cm.plasma(np.linspace(0,1,len(Ls)))
for idx,L in enumerate(Ls):
    xs, ys, ses = by_L[L]
    mask = ~np.isnan(ys)
    if np.sum(mask) < 2:
        print(f'  Skipping L={L}: insufficient data')
        continue
    ax.plot(xs[mask], ys[mask], '-o', label=f'L={L}', color=colors[idx])
    ax.fill_between(xs[mask], ys[mask]-ses[mask], ys[mask]+ses[mask], alpha=0.15, color=colors[idx])

ax.set_xlabel('Coupling ratio (Γ)', fontsize=12)
ax.set_ylabel('Binder cumulant U', fontsize=12)
ax.set_title('Binder U vs Γ — L=48,32,24,16 (15 ratios)', fontsize=13)
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plot_path = bench / 'binder_range_Ls.png'
plt.savefig(plot_path, dpi=150)
print(f'Plot saved to {plot_path}')

# Print short table to stdout
for L in Ls:
    xs, ys, ses = by_L[L]
    valid = ~np.isnan(ys)
    print(f'\nL={L}: {np.sum(valid)}/15 ratios present')
    for x,y,se in zip(xs, ys, ses):
        if np.isnan(y):
            print(f'  γ={x:.6f}: MISSING')
        else:
            print(f'  γ={x:.6f}: U={y:.5f} ± {se:.5f}')
