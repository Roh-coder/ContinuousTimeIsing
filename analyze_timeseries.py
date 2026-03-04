#!/usr/bin/env python3
"""
Analyze binder_Ls_timeseries.csv: compute statistics, autocorr times, and precision goals.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import colorsys

csv.field_size_limit(sys.maxsize)

def tau_int(x, max_lag=None):
    """Integrated autocorrelation time via FFT."""
    n = len(x)
    if n < 2:
        return 1.0
    x = x - np.mean(x)
    if max_lag is None:
        max_lag = min(n - 1, int(n // 2))
    
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    if acf[0] == 0:
        return 1.0
    acf /= acf[0]
    
    tau = 0.5
    for t in range(1, min(len(acf), n // 2)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return max(tau, 0.5)

bench = Path('bench')
ts_file = bench / 'binder_Ls_timeseries.csv'

if not ts_file.exists():
    print(f'Error: {ts_file} not found')
    sys.exit(1)

# Load timeseries data
print(f'Reading {ts_file}...')
data_by_Lratio = defaultdict(list)
total_rows = 0

with open(ts_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_rows += 1
        try:
            L = int(row['L'])
            ratio = float(row['ratio'])
            U = float(row['U'])
            if 0 <= U <= 1:  # Valid U values
                key = (L, round(ratio, 6))
                data_by_Lratio[key].append(U)
        except Exception:
            pass

print(f'  Rows: {total_rows}, (L,γ) pairs: {len(data_by_Lratio)}')

# Compute summary statistics
summary = []
for (L, ratio), U_list in sorted(data_by_Lratio.items()):
    if len(U_list) < 10:
        continue
    
    U_arr = np.array(U_list)
    mean_U = float(np.mean(U_arr))
    var_U = float(np.var(U_arr, ddof=1))
    se_U = float(np.std(U_arr, ddof=1) / np.sqrt(len(U_arr)))
    
    # Autocorrelation time
    tau = tau_int(U_arr)
    
    # Configs needed for 1% precision (δU = 0.01)
    delta_target = 0.01
    # se_U = sqrt(var_U * tau) / sqrt(n_cfg)
    # n_cfg = var_U * tau / (delta_target^2)
    if var_U > 0:
        n_configs_1pct = int(np.ceil(var_U * tau / (delta_target**2)))
        n_autocorr_1pct = int(np.ceil(n_configs_1pct / tau))
    else:
        n_configs_1pct = 1
        n_autocorr_1pct = 1
    
    # Configs needed for 0.1% precision
    delta_target_001 = 0.001
    if var_U > 0:
        n_configs_001pct = int(np.ceil(var_U * tau / (delta_target_001**2)))
        n_autocorr_001pct = int(np.ceil(n_configs_001pct / tau))
    else:
        n_configs_001pct = 1
        n_autocorr_001pct = 1
    
    summary.append({
        'L': L, 'ratio': ratio, 'n': len(U_list),
        'mean_U': mean_U, 'var_U': var_U, 'se_U': se_U,
        'tau_int': tau,
        'n_cfg_1pct': n_configs_1pct,
        'n_cfg_001pct': n_configs_001pct,
        'n_tau_1pct': n_autocorr_1pct,
        'n_tau_001pct': n_autocorr_001pct
    })

# Write summary CSV
summary_file = bench / 'summary_binder_Ls.csv'
with open(summary_file, 'w') as f:
    f.write('L,ratio,n_measured,mean_U,var_U,se_U,tau_int,n_configs_1pct,n_autocorr_1pct,n_configs_001pct,n_autocorr_001pct\n')
    for row in summary:
        f.write(f"{row['L']},{row['ratio']:.10f},{row['n']},"
                f"{row['mean_U']:.8f},{row['var_U']:.8f},{row['se_U']:.8f},"
                f"{row['tau_int']:.2f},{row['n_cfg_1pct']},{row['n_tau_1pct']},"
                f"{row['n_cfg_001pct']},{row['n_tau_001pct']}\n")

print(f'Summary: {summary_file} ({len(summary)} points)')

# Group by L for plotting
by_L = defaultdict(list)
for row in summary:
    by_L[row['L']].append((row['ratio'], row['mean_U'], row['se_U']))

# Plot
fig, ax = plt.subplots(figsize=(12, 7))
hues = np.linspace(0.0, 1.0, len(by_L), endpoint=False)
colors = [colorsys.hls_to_rgb(h, 0.26, 0.75) for h in hues]

for idx, L in enumerate(sorted(by_L.keys())):
    data = sorted(by_L[L], key=lambda x: x[0])
    if len(data) < 2:
        continue
    
    ratios, Us, ses = zip(*data)
    ratios, Us, ses = np.array(ratios), np.array(Us), np.array(ses)
    
    # Add uncertainty ribbons FIRST (so they appear behind)
    ax.fill_between(ratios, Us - ses, Us + ses, alpha=0.45, color=colors[idx], 
                    label=f'L={L} (±1σ)')
    
    # Add error bars at each point with thin cross markers
    ax.errorbar(ratios, Us, yerr=ses, fmt='x', color=colors[idx], 
                elinewidth=2, capsize=5, markersize=9, markeredgewidth=1.5,
                alpha=0.85, capthick=2)
    
    # Then plot thin line on top connecting points
    ax.plot(ratios, Us, '-', linewidth=1.3, color=colors[idx], alpha=0.7)

ax.set_xlabel('Coupling ratio (Γ)', fontsize=13, fontweight='bold')
ax.set_ylabel('Binder cumulant U', fontsize=13, fontweight='bold')
ax.set_title('Binder Crossing Diagram with Uncertainties: L=48,32,24,16 (15 log-spaced ratios)', 
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='best', framealpha=0.95, title='System size (±1σ errors)')
ax.set_xscale('log')
ax.grid(True, alpha=0.3, linestyle='--', which='both')
ax.axvline(x=1.0, color='red', linestyle=':', linewidth=2, alpha=0.4)

plt.tight_layout()
plot_file = bench / 'binder_crossing_final.png'
plt.savefig(plot_file, dpi=150, bbox_inches='tight')
print(f'Plot: {plot_file}')

# Print summary table
print('\n' + '='*130)
print('DETAILED SUMMARY TABLE (with autocorrelation time units)')
print('='*130)
for L in sorted(by_L.keys()):
    points = sorted(by_L[L])
    summary_for_L = [s for s in summary if s['L'] == L]
    print(f'\nL={L:2d} ({len(points)} ratios):')
    print('-' * 130)
    print(f'{"Ratio":>8} | {"U":>7} | {"±δU":>7} | {"Var":>8} | {"τ_int":>6} | {"Cfgs(1%)":>10} | {"τ(1%)":>8} | {"Cfgs(0.1%)":>12} | {"τ(0.1%)":>10}')
    print('-' * 130)
    for s in summary_for_L:
        ratio = s['ratio']
        mean_U = s['mean_U']
        se_U = s['se_U']
        var_U = s['var_U']
        tau = s['tau_int']
        n_1pct = s['n_cfg_1pct']
        n_tau_1pct = s['n_tau_1pct']
        n_001pct = s['n_cfg_001pct']
        n_tau_001pct = s['n_tau_001pct']
        print(f'{ratio:8.6f} | {mean_U:7.4f} | {se_U:7.5f} | {var_U:8.5f} | {tau:6.2f} | {n_1pct:10d} | {n_tau_1pct:8d} | {n_001pct:12d} | {n_tau_001pct:10d}')

print('\n' + '='*130)
