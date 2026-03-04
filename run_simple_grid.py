#!/usr/bin/env python3
"""
Simple grid: L=48,32,24,16 × 15 log-spaced γ in (0.9, 10/9), 500 configs each.
Fixed: no loops, explicit job list, clear tracking.
"""
import subprocess
import numpy as np
from pathlib import Path
import sys

L_vals = [48, 32, 24, 16]
n_ratios = 15
n_configs = 500

ratios = np.logspace(np.log10(0.9), np.log10(10.0/9.0), n_ratios)

bench_dir = Path('bench')
bench_dir.mkdir(exist_ok=True)

con_exe = Path('build/ConTimeIsing')
if not con_exe.exists():
    print('Error: build/ConTimeIsing not found')
    sys.exit(1)

out_csv = bench_dir / 'binder_Ls_timeseries.csv'
with open(out_csv, 'w') as f:
    f.write('L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n')

# Build job list explicitly
jobs = []
for L in L_vals:
    for idx, ratio in enumerate(ratios):
        jobs.append((L, ratio, idx))

print(f'Running {len(jobs)} jobs: 4 Ls × 15 ratios × 500 configs')

for job_idx, (L, ratio, r_idx) in enumerate(jobs, 1):
    tmp_file = bench_dir / f'tmp_L{L}_r{ratio:.8f}.csv'
    
    # Clean up any stale tmp file
    if tmp_file.exists():
        tmp_file.unlink()
    
    cmd = [str(con_exe),
           '--layers', str(L),
           '--ratio_min', f'{ratio:.10f}',
           '--ratio_max', f'{ratio:.10f}',
           '--n_configs', str(n_configs),
           '--burn_in', '250',
           '--sweeps_between', '5',
           '--seed', str(99000 + L*100 + r_idx),
           '--ts_file', str(tmp_file)]
    
    print(f'[{job_idx:2d}/60] L={L:2d} γ={ratio:.6f}', end=' ', flush=True)
    
    try:
        result = subprocess.run(cmd, check=False, timeout=300, capture_output=True)
        if result.returncode != 0:
            print('✗ ret != 0')
            if tmp_file.exists():
                tmp_file.unlink()
            continue
    except subprocess.TimeoutExpired:
        print('✗ timeout')
        if tmp_file.exists():
            tmp_file.unlink()
        continue
    except Exception as e:
        print(f'✗ {e}')
        if tmp_file.exists():
            tmp_file.unlink()
        continue
    
    # Append to output CSV
    if tmp_file.exists():
        try:
            lines = 0
            with open(tmp_file) as tf, open(out_csv, 'a') as of:
                for line in tf:
                    if line.strip() and not line.startswith('L,'):
                        of.write(line)
                        lines += 1
            tmp_file.unlink()
            print(f'✓ {lines} rows')
        except Exception as e:
            print(f'✗ append failed: {e}')
    else:
        print('✗ tmp file missing')

print(f'\nDone: {out_csv}')
sys.exit(0)
