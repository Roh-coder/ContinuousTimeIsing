#!/usr/bin/env python3
import csv, numpy as np
from pathlib import Path

def tau_int(x):
    n = len(x)
    if n < 2: return 0.0
    x = x - np.mean(x)
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for t in range(1, min(len(acf), n//2)):
        if acf[t] <= 0: break
        tau += acf[t]
    return tau

results = []
for ts_file in sorted(Path('bench').glob('ts_L*.csv')):
    rows = []
    with open(ts_file) as f:
        csv_rdr = csv.reader(f)
        next(csv_rdr)
        for row in csv_rdr:
            if len(row) > 6:
                rows.append({'L': int(row[0]), 'ratio': float(row[1]), 'U': float(row[6])})
    if not rows: 
        continue
    groups = {}
    for r in rows:
        key = r['ratio']
        if key not in groups: 
            groups[key] = []
        groups[key].append(r['U'])
    L = rows[0]['L']
    for ratio in sorted(groups.keys()):
        Uarr = np.array(groups[ratio])
        mean_U = np.mean(Uarr)
        var_U = np.var(Uarr, ddof=1)
        t = tau_int(Uarr)
        ac_u = var_U / (0.001**2)
        results.append((L, ratio, len(Uarr), mean_U, var_U, t, ac_u))

with open('bench/performance_summary.csv', 'w') as f:
    f.write("L,ratio,n_samples,mean_U,var_U,tau_int_configs,autocorr_times_for_eps_001\n")
    for r in results:
        f.write(f"{int(r[0])},{r[1]:.1f},{int(r[2])},{r[3]:.6f},{r[4]:.6f},{r[5]:.2f},{r[6]:.1f}\n")

print("Summary written to bench/performance_summary.csv")
