#!/usr/bin/env python3
import csv, numpy as np
from pathlib import Path
from collections import defaultdict

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

# Open output file
out = open('bench/summary_per_point.csv', 'w')
out.write('L,ratio,n_measured,mean_U,var_U,tau_int (configs),N_required (configs),ac_units (N/(2*tau))\n')

# Process each ts file
eps = 0.001
for ts_path in sorted(Path('bench').glob('ts_L*.csv')):
    rows = []
    with open(ts_path) as f:
        csv_rdr = csv.reader(f)
        next(csv_rdr)  # skip header
        for row in csv_rdr:
            if len(row) > 6:
                rows.append({'L': int(row[0]), 'ratio': float(row[1]), 'U': float(row[6])})
    
    if not rows: continue
    groups = defaultdict(list)
    for r in rows:
        groups[r['ratio']].append(r['U'])
    
    L = rows[0]['L']
    for ratio in sorted(groups.keys()):
        Uarr = np.array(groups[ratio])
        mean_U = float(np.mean(Uarr))
        var_U = float(np.var(Uarr, ddof=1))
        t = tau_int(Uarr)
        N_req = (2.0 * t * var_U) / (eps**2)
        ac_u = var_U / (eps**2)
        out.write(f'{L},{ratio:.2f},{len(Uarr)},{mean_U:.10g},{var_U:.10g},{t:.6g},{N_req:.6g},{ac_u:.6g}\n')

out.close()
print('Done')
