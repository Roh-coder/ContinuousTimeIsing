#!/usr/bin/env python3
"""Analyze existing timeseries files and produce summary."""
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np

def integrated_autocorrelation_time(x, max_lag=None):
    n = len(x)
    if n < 2:
        return 0.0
    x = x - np.mean(x)
    if max_lag is None:
        max_lag = min(n - 1, int(n / 2))
    
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    acf /= acf[0]
    acf = acf[: max_lag + 1]
    
    tau = 0.5
    for t in range(1, len(acf)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return float(tau)

def analyze(U, eps=1e-3):
    n = len(U)
    if n < 2:
        return None
    mean = float(np.mean(U))
    var = float(np.var(U, ddof=1))
    tau = integrated_autocorrelation_time(U)
    N_required = (2.0 * tau * var) / (eps * eps) if eps > 0 else float('inf')
    ac_units = var / (eps * eps)
    return {
        'n': n,
        'mean': mean,
        'var': var,
        'tau': tau,
        'N_required': N_required,
        'ac_units': ac_units,
    }

# Process all ts files
results = []
for ts_file in sorted(Path('bench').glob('ts_L*.csv')):
    print(f"\nProcessing {ts_file.name}...")
    rows = []
    with open(ts_file) as f:
        rdr = csv.reader(f)
        try:
            header = next(rdr)
        except StopIteration:
            continue
        for row in rdr:
            if not row:
                continue
            try:
                L = int(row[0])
                ratio = float(row[1])
                U = float(row[6])
                rows.append({'L': L, 'ratio': ratio, 'U': U})
            except:
                pass
    
    # Group by ratio
    groups = defaultdict(list)
    for r in rows:
        groups[r['ratio']].append(r['U'])
    
    L = rows[0]['L'] if rows else 0
    for ratio, Ulist in sorted(groups.items()):
        Uarr = np.array(Ulist, dtype=float)
        info = analyze(Uarr, eps=0.001)
        if info:
            print(f"  L={L} ratio={ratio:6.2f}: n={info['n']:4d} mean_U={info['mean']:8.5f} var={info['var']:10.6f} tau={info['tau']:8.2f} ac_units={info['ac_units']:10.1f}")
            results.append({
                'L': L,
                'ratio': ratio,
                'n_measured': int(info['n']),
                'mean_U': info['mean'],
                'var_U': info['var'],
                'tau_int (configs)': info['tau'],
                'N_required (configs)': info['N_required'],
                'ac_units (N/(2*tau))': info['ac_units'],
            })

# Save summary
keys = ['L','ratio','n_measured','mean_U','var_U','tau_int (configs)','N_required (configs)','ac_units (N/(2*tau))']
with open('bench/summary_per_point.csv', 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=keys)
    wr.writeheader()
    for r in sorted(results, key=lambda x: (x['L'], x['ratio'])):
        wr.writerow(r)

print(f"\n\nSummary saved to bench/summary_per_point.csv")
print(f"Total points: {len(results)}")
