#!/usr/bin/env python3
"""
High-resolution benchmark with tight ratios around criticality.
- Lattice sizes: L = 4, 8, 12, 16, 20, 24, 28, 32
- Ratios: dense scan from 0.90 to 1.11
- Jobs scheduled: highest cost first (largest L, highest ratio)
- Output: Single consolidated timeseries CSV
"""
import csv
import numpy as np
import subprocess
import sys
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

def parse_ts_file(path):
    """Parse timeseries CSV."""
    rows = []
    with open(path) as f:
        reader = csv.reader(f)
        try:
            next(reader)  # skip header
        except StopIteration:
            return rows
        for row in reader:
            if row:
                try:
                    rows.append({
                        'L': int(row[0]),
                        'ratio': float(row[1]),
                        'U': float(row[6])
                    })
                except (IndexError, ValueError):
                    pass
    return rows

# Configuration
L_min, L_max = 4, 32
ratio_min, ratio_max = 0.90, 1.11
n_configs = 1000  # High accuracy
burn_in = 300
sweeps_between = 5
seed_base = 12345

# Generate lattice sizes: 4, 8, 12, 16, 20, 24, 28, 32
Ls = list(range(L_min, L_max + 1, 4))
print(f"Lattice sizes: {Ls}")

# Generate dense ratio grid around critical point
n_ratios = 22  # Dense scan
ratios = np.linspace(ratio_min, ratio_max, n_ratios)
print(f"Ratio points: {len(ratios)}, from {ratios[0]:.4f} to {ratios[-1]:.4f}")

# Generate all jobs
jobs = []
for L in Ls:
    for ratio in ratios:
        # Estimated cost: proportional to L and distance from criticality
        # Higher L and ratios near 1.0 are more expensive
        criticality_cost = 1.0 / (0.1 + abs(ratio - 1.0))  # explodes near 1.0
        cost = L * criticality_cost
        jobs.append({'L': L, 'ratio': ratio, 'cost': cost})

# Sort by cost (descending) - expensive jobs first
jobs.sort(key=lambda x: x['cost'], reverse=True)

print(f"\nTotal jobs: {len(jobs)}")
print(f"Most expensive: L={jobs[0]['L']}, ratio={jobs[0]['ratio']:.4f}, cost={jobs[0]['cost']:.1f}")
print(f"Least expensive: L={jobs[-1]['L']}, ratio={jobs[-1]['ratio']:.4f}, cost={jobs[-1]['cost']:.1f}")

# Run all jobs, consolidate output
print("\nRunning high-resolution benchmark...\n")

ts_output_file = Path("bench") / "ts_high_res.csv"
ts_output_file.parent.mkdir(exist_ok=True)

# Write header
with open(ts_output_file, 'w') as f:
    f.write("L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n")

total_jobs = len(jobs)
all_data = defaultdict(list)

for job_idx, job in enumerate(jobs):
    L = job['L']
    ratio = job['ratio']
    cost = job['cost']
    
    # Create temporary timeseries file
    tmp_ts_file = Path("bench") / f"tmp_ts_L{L}_r{ratio:.4f}.csv"
    
    # Build command
    cmd = [
        "./build/ConTimeIsing",
        "--layers", str(L),
        "--n_configs", str(n_configs),
        "--sweeps_between", str(sweeps_between),
        "--burn_in", str(burn_in),
        "--seed", str(seed_base + job_idx),
        "--ratio_min", str(ratio),
        "--ratio_max", str(ratio),
        "--ratio_n", "1",
        "--ts_file", str(tmp_ts_file),
    ]
    
    print(f"[{job_idx+1}/{total_jobs}] L={L:2d}, ratio={ratio:.4f}, cost={cost:.1f}", end=" ... ")
    sys.stdout.flush()
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"ERROR: {result.stderr[:100]}")
            continue
        
        # Parse and append to main file
        rows = parse_ts_file(tmp_ts_file)
        with open(ts_output_file, 'a') as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow([
                    row['L'], row['ratio'], 
                    # Fill in dummy cfg and other fields from the original row
                    0, 0, 0, 0, row['U'], 0
                ])
        
        # Track for summary
        for row in rows:
            all_data[(L, ratio)].append(row['U'])
        
        # Clean up
        tmp_ts_file.unlink(missing_ok=True)
        print(f"✓ {len(rows)} configs")
        
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT")
        tmp_ts_file.unlink(missing_ok=True)
    except Exception as e:
        print(f"FAIL: {e}")
        tmp_ts_file.unlink(missing_ok=True)

# Generate summary
print("\n" + "="*80)
print("SUMMARY (High-Res Scan: L=4-32, ratio=0.9-1.11)")
print("="*80 + "\n")

summary_data = []
for L in Ls:
    print(f"L={L}:")
    for ratio in ratios:
        if (L, ratio) in all_data:
            U_vals = np.array(all_data[(L, ratio)])
            mean_U = np.mean(U_vals)
            std_U = np.std(U_vals, ddof=1) if len(U_vals) > 1 else 0
            tau = tau_int(U_vals)
            var_U = np.var(U_vals, ddof=1) if len(U_vals) > 1 else 0
            ac_units = var_U / (0.001**2)
            
            summary_data.append([L, ratio, len(U_vals), mean_U, var_U, tau, ac_units])
            print(f"  γ={ratio:.4f}: U={mean_U:.5f}±{std_U:.5f}, τ={tau:.2f}, ac_units={ac_units:.0f}")

# Write summary CSV
summary_file = Path("bench") / "summary_high_res.csv"
with open(summary_file, 'w') as f:
    f.write("L,ratio,n_measured,mean_U,var_U,tau_int,ac_units\n")
    writer = csv.writer(f)
    writer.writerows(summary_data)

print(f"\n✓ Timeseries saved to: {ts_output_file}")
print(f"✓ Summary saved to: {summary_file}")
print(f"✓ Total data points: {len(summary_data)}")
