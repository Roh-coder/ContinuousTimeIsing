#!/usr/bin/env python3
"""
Generate high-res summary and visualization from ts_high_res.csv
Runs after benchmark completion.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

def tau_int(x, max_lag=None):
    """Integrated autocorrelation time via FFT."""
    n = len(x)
    if n < 2:
        return 0.0
    x = x - np.mean(x)
    if max_lag is None:
        max_lag = min(n - 1, int(n / 2))
    
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    acf /= acf[0]
    
    tau = 0.5
    for t in range(1, min(len(acf), n // 2)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return tau

# Read timeseries
ts_path = Path("bench/ts_high_res.csv")
if not ts_path.exists():
    print(f"File not found: {ts_path}")
    exit(1)

data_by_L_ratio = defaultdict(list)
import sys
csv.field_size_limit(sys.maxsize)
with open(ts_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            L = int(row['L'])
            ratio = float(row['ratio'])
            U = float(row['U'])
            data_by_L_ratio[(L, ratio)].append(U)
        except (ValueError, KeyError):
            pass

print(f"Processed {len(data_by_L_ratio)} (L, ratio) points")

# Compute statistics
summary = []
for (L, ratio), U_list in sorted(data_by_L_ratio.items()):
    if len(U_list) < 10:
        continue
    
    U_arr = np.array(U_list)
    mean_U = float(np.mean(U_arr))
    var_U = float(np.var(U_arr, ddof=1))
    se_U = np.sqrt(var_U / len(U_arr))
    tau = tau_int(U_arr)
    ac_units = var_U / (0.001**2)
    
    summary.append({
        'L': L, 'ratio': ratio,
        'mean_U': mean_U, 'var_U': var_U, 'se_U': se_U,
        'tau': tau, 'ac_units': ac_units, 'n': len(U_list)
    })

# Write summary CSV
summary_path = Path("bench/summary_high_res.csv")
with open(summary_path, 'w') as f:
    f.write("L,ratio,n_measured,mean_U,var_U,se_U,tau_int,ac_units\n")
    for row in summary:
        f.write(f"{row['L']},{row['ratio']},{row['n']},{row['mean_U']:.6f},"
                f"{row['var_U']:.6f},{row['se_U']:.6f},"
                f"{row['tau']:.2f},{row['ac_units']:.0f}\n")

print(f"Summary written to {summary_path}")

# Group by L for plotting
by_L = defaultdict(list)
for row in summary:
    by_L[row['L']].append((row['ratio'], row['mean_U'], row['se_U']))

# Plot with error ribbons
fig, ax = plt.subplots(figsize=(14, 9))
colors = plt.cm.viridis(np.linspace(0, 1, len(by_L)))

for idx, L in enumerate(sorted(by_L.keys())):
    data = sorted(by_L[L], key=lambda x: x[0])
    if len(data) < 2:
        continue
    
    ratios, Us, ses = zip(*data)
    ratios, Us, ses = np.array(ratios), np.array(Us), np.array(ses)
    
    # Line + markers
    ax.plot(ratios, Us, 'o-', linewidth=2.5, markersize=6,
            label=f'L={L}', color=colors[idx])
    
    # Error ribbons
    upper = Us + ses
    lower = Us - ses
    ax.fill_between(ratios, lower, upper, alpha=0.15, color=colors[idx])

ax.set_xlabel('Coupling Ratio (Γ)', fontsize=14, fontweight='bold')
ax.set_ylabel('Binder Cumulant U', fontsize=14, fontweight='bold')
ax.set_title('High-Resolution Binder Crossing (L=4-32, Dense Ratio Scan)',
             fontsize=15, fontweight='bold')
ax.legend(fontsize=11, loc='best', framealpha=0.95, ncol=2)
ax.grid(True, alpha=0.3, linestyle='--')
ax.axvline(x=1.0, color='red', linestyle=':', linewidth=2, alpha=0.5, label='γ=1.0')

plt.tight_layout()
plot_path = Path("bench/binder_high_res.png")
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {plot_path}")

# Print summary table
print("\n" + "="*80)
print("HIGH-RES BENCHMARK SUMMARY")
print("="*80 + "\n")
for L in sorted(by_L.keys()):
    print(f"L={L}: {len(by_L[L])} ratio points")
    for ratio, U, se in sorted(by_L[L]):
        print(f"  γ={ratio:.4f}: U={U:.5f}±{se:.5f}")
