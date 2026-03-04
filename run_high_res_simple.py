#!/usr/bin/env python3
"""
High-resolution benchmark: L=4-32, dense ratio scan 9/10 to 10/9
500 configs each, jobs scheduled expensive-first (L descending, ratio descending)
Output: Single consolidated timeseries CSV
"""
import subprocess
import sys
from pathlib import Path

# Configuration
L_vals = [32, 28, 24, 20, 16, 12, 8, 4]  # Descending (expensive first)
n_ratios = 21  # Dense scan
ratio_min = 0.9           # 9/10
ratio_max = 10.0/9.0      # 10/9 ≈ 1.1111
n_configs = 500
burn_in = 250
sweeps_between = 5
seed_base = 54321

# Generate ratio grid (high to low for scheduling)
import numpy as np
ratios = np.linspace(ratio_max, ratio_min, n_ratios)  # High to low

print(f"Configuration:")
print(f"  Lattice sizes: {L_vals} (high to low)")
print(f"  Ratios: {n_ratios} points from {ratio_min:.4f} to {ratio_max:.4f}")
print(f"  Per-point: {n_configs} configs, burn_in={burn_in}, sweeps_between={sweeps_between}")
print(f"  Total jobs: {len(L_vals) * n_ratios}")
print()

# Create output file with header
ts_file = Path("bench/ts_high_res.csv")
ts_file.parent.mkdir(exist_ok=True)
with open(ts_file, 'w') as f:
    f.write("L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n")

# Run jobs: L descending, then ratio descending (expensive first)
total_jobs = len(L_vals) * len(ratios)
job_idx = 0

for L in L_vals:
    for ratio in ratios:
        job_idx += 1
        
        # Create temp file for this job
        tmp_file = Path("bench") / f"tmp_L{L}_r{ratio:.6f}.csv"
        
        # Build command
        cmd = [
            "./build/ConTimeIsing",
            "--layers", str(L),
            "--ratio_min", f"{ratio:.10f}",
            "--ratio_max", f"{ratio:.10f}",
            "--ratio_n", "1",
            "--n_configs", str(n_configs),
            "--sweeps_between", str(sweeps_between),
            "--burn_in", str(burn_in),
            "--seed", str(seed_base + job_idx),
            "--ts_file", str(tmp_file),
            "--print_every_burn", "0",
            "--print_every_cfg", "0",
        ]
        
        print(f"[{job_idx:4d}/{total_jobs}] L={L:2d}, γ={ratio:.6f} ...", end=" ", flush=True)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                print(f"ERROR")
                continue
            
            # Append to main file
            if tmp_file.exists():
                with open(tmp_file) as f_in:
                    lines = f_in.readlines()
                    # Skip header, append data
                    with open(ts_file, 'a') as f_out:
                        f_out.writelines(lines[1:])
                tmp_file.unlink()
                print(f"✓")
            else:
                print(f"NO OUTPUT")
        except subprocess.TimeoutExpired:
            print(f"TIMEOUT")
            tmp_file.unlink(missing_ok=True)
        except Exception as e:
            print(f"FAIL: {e}")
            tmp_file.unlink(missing_ok=True)

print(f"\n✓ All data saved to: {ts_file}")
print(f"  Run `python3 analyze_high_res.py` to generate summary and plots")
