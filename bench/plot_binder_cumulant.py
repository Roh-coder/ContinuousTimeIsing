#!/usr/bin/env python3
"""
Generate a Binder cumulant crossing plot from summary CSV data.
"""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Read the CSV file
csv_file = '/workspaces/ContinuousTimeIsing/bench/summary_per_point.csv'
data_by_L = defaultdict(lambda: {'ratio': [], 'mean_U': []})

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        L = int(row['L'])
        ratio = float(row['ratio'])
        mean_U = float(row['mean_U'])
        data_by_L[L]['ratio'].append(ratio)
        data_by_L[L]['mean_U'].append(mean_U)

# Sort the data by ratio for each L
for L in data_by_L:
    # Sort by ratio, keeping corresponding mean_U values in sync
    sorted_pairs = sorted(zip(data_by_L[L]['ratio'], data_by_L[L]['mean_U']))
    data_by_L[L]['ratio'] = [p[0] for p in sorted_pairs]
    data_by_L[L]['mean_U'] = [p[1] for p in sorted_pairs]

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each L value
for L in sorted(data_by_L.keys()):
    ax.plot(data_by_L[L]['ratio'], data_by_L[L]['mean_U'], marker='o', label=f'L={L}')

# Configure plot appearance
ax.set_xlabel('Coupling Ratio (Gamma)', fontsize=12)
ax.set_ylabel('Binder Cumulant U', fontsize=12)
ax.set_title('Binder Cumulant Crossing Plot', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)

# Save the plot
output_file = '/workspaces/ContinuousTimeIsing/bench/binder_vs_ratio.png'
plt.savefig(output_file, dpi=120, bbox_inches='tight')
print(f"Plot saved to: {output_file}")

plt.close()
