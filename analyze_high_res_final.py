#!/usr/bin/env python3
"""
Analyze high-res benchmark data: generate summary CSV and Binder crossing plot
"""
import csv
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

# Increase field size limit
csv.field_size_limit(sys.maxsize)

def tau_int(x, max_lag=None):
    """Integrated autocorrelation time via FFT."""
    n = len(x)
    if n < 2:
        return 0.0
    x = x - np.mean(x)
    if max_lag is None:
        max_lag = min(n - 1, int(n // 2))
    
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    if acf[0] == 0:
        return 0.0
    acf /= acf[0]
    
    tau = 0.5
    for t in range(1, min(len(acf), n // 2)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return tau

print("Reading timeseries data...")
ts_path = Path("bench/ts_high_res.csv")
if not ts_path.exists():
    print(f"ERROR: {ts_path} not found")
    sys.exit(1)

data_by_L_ratio = defaultdict(list)
total_rows = 0

with open(ts_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_rows += 1
        try:
            L = int(row['L'])
            ratio = float(row['ratio'])
            U = float(row['U'])
            if U > 0:  # Skip uninitialized/dummy values
                data_by_L_ratio[(L, ratio)].append(U)
        except (ValueError, KeyError):
            pass

print(f"  Rows processed: {total_rows}")
print(f"  Unique (L, γ) points: {len(data_by_L_ratio)}")

# Compute statistics
summary = []
for (L, ratio), U_list in sorted(data_by_L_ratio.items()):
    if len(U_list) < 5:
        continue
    
    U_arr = np.array(U_list)
    mean_U = float(np.mean(U_arr))
    var_U = float(np.var(U_arr, ddof=1))
    se_U = np.sqrt(var_U / len(U_arr))
    tau = tau_int(U_arr)
    ac_units = var_U / (0.001**2) if var_U > 0 else 0
    
    summary.append({
        'L': L, 'ratio': ratio, 'n': len(U_list),
        'mean_U': mean_U, 'var_U': var_U, 'se_U': se_U,
        'tau': tau, 'ac_units': ac_units
    })

# Write summary CSV
summary_path = Path("bench/summary_high_res_final.csv")
with open(summary_path, 'w') as f:
    f.write("L,ratio,n_measured,mean_U,var_U,se_U,tau_int,ac_units\n")
    for row in summary:
        f.write(f"{row['L']},{row['ratio']:.10f},{row['n']},{row['mean_U']:.6f},"
                f"{row['var_U']:.6f},{row['se_U']:.6f},"
                f"{row['tau']:.2f},{row['ac_units']:.0f}\n")

print(f"Summary written: {summary_path} ({len(summary)} points)")

# Group by L for plotting
by_L = defaultdict(list)
for row in summary:
    by_L[row['L']].append((row['ratio'], row['mean_U'], row['se_U']))

if not by_L:
    print("No data to plot")
    sys.exit(0)

# Create plot
print("Creating plot...")
fig, ax = plt.subplots(figsize=(14, 9))
colors = plt.cm.viridis(np.linspace(0, 1, len(by_L)))

for idx, L in enumerate(sorted(by_L.keys())):
    data = sorted(by_L[L], key=lambda x: x[0])
    if len(data) < 2:
        continue
    
    ratios, Us, ses = zip(*data)
    ratios, Us, ses = np.array(ratios), np.array(Us), np.array(ses)
    
    ax.plot(ratios, Us, 'o-', linewidth=2.5, markersize=7,
            label=f'L={L}', color=colors[idx])
    
    upper = Us + ses
    lower = Us - ses
    ax.fill_between(ratios, lower, upper, alpha=0.15, color=colors[idx])

ax.set_xlabel('Coupling Ratio (Γ)', fontsize=14, fontweight='bold')
ax.set_ylabel('Binder Cumulant U', fontsize=14, fontweight='bold')
ax.set_title('High-Resolution Binder Crossing (L=4-32, γ=0.9-10/9, 500 configs each)',
             fontsize=15, fontweight='bold')
ax.legend(fontsize=11, loc='best', framealpha=0.95, ncol=2)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axvline(x=1.0, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax.text(1.0, ax.get_ylim()[1]*0.95, 'γ=1.0', fontsize=10, ha='right', color='red')

plt.tight_layout()
plot_path = Path("bench/binder_high_res_final.png")
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
print(f"Plot saved: {plot_path}")

# Summary table
print("\n" + "="*100)
print("HIGH-RESOLUTION BENCHMARK SUMMARY")
print("="*100 + "\n")
for L in sorted(by_L.keys()):
    points = sorted(by_L[L])
    print(f"L={L:2d} ({len(points)} ratio points):")
    for ratio, U, se in points:
        print(f"  γ={ratio:.6f}: U={U:.5f}±{se:.5f}")
