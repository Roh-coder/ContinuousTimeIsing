#!/usr/bin/env python3
"""Generate Binder cumulant crossing plot."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Read summary
Ls = {}
with open('summary_per_point.csv') as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        L = int(row['L'])
        ratio = float(row['ratio'])
        mean_U = float(row['mean_U'])
        if L not in Ls:
            Ls[L] = []
        Ls[L].append((ratio, mean_U))

# Plot
fig, ax = plt.subplots(figsize=(10, 7))
for L in sorted(Ls.keys()):
    xs, ys = zip(*sorted(Ls[L]))
    ax.plot(xs, ys, 'o-', linewidth=2.5, markersize=10, label=f'L={L}')

ax.set_xlabel('Coupling ratio (Gamma)', fontsize=13, fontweight='bold')
ax.set_ylabel('Binder Cumulant U', fontsize=13, fontweight='bold')
ax.set_title('Binder Cumulant Crossing\n(Continuous-Time Ising Model)', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.4, linestyle='--')
ax.set_ylim([0.55, 0.75])

plt.tight_layout()
plt.savefig('binder_vs_ratio.png', dpi=120, bbox_inches='tight')
print("✓ Plot saved to binder_vs_ratio.png")
