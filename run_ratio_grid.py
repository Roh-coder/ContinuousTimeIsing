#!/usr/bin/env python3
"""
Targeted runner for L=[48,32,24,16], 15 ratios (10/9 to 0.9), 1000 configs per point.
Skips ratios already in per-L CSVs, writes to bench/ts_L<N>_high_res.csv
"""
import subprocess
import sys
import time
from pathlib import Path
import csv
import numpy as np

L_vals = [48, 32, 24, 16]
ratio_min = 0.9
ratio_max = 10.0 / 9.0
n_ratios = 15
n_configs = 1000
burn_in = 250
sweeps_between = 5
seed_base = 54321

ratios = list(reversed([ratio_min + i*(ratio_max-ratio_min)/(n_ratios-1) for i in range(n_ratios)]))

bench_dir = Path('bench')
bench_dir.mkdir(exist_ok=True)

con_exe = Path('build/ConTimeIsing')
if not con_exe.exists():
    print('Error: build/ConTimeIsing not found. Build the project first.'); sys.exit(1)

for L in L_vals:
    out_file = bench_dir / f'ts_L{L}_high_res.csv'
    existing_ratios = set()
    
    if out_file.exists():
        with open(out_file) as f:
            reader = csv.DictReader(f)
            for r in reader:
                try:
                    existing_ratios.add(round(float(r['ratio']), 6))
                except Exception:
                    pass
    else:
        with open(out_file, 'w') as f:
            f.write('L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n')

    print(f'\nRunning L={L} ({len(ratios)} ratios × {n_configs} configs)')
    print(f'  Existing: {len(existing_ratios)} ratios, skipping those')

    for idx, ratio in enumerate(ratios):
        ratio_rounded = round(ratio, 6)
        if ratio_rounded in existing_ratios:
            print(f'  [{idx+1:2d}/{len(ratios)}] γ={ratio:.6f}: SKIP (exists)')
            continue

        tmp_file = bench_dir / f'tmp_ts_L{L}_r{ratio:.6f}.csv'
        cmd = [str(con_exe),
               '--layers', str(L),
               '--ratio_min', f'{ratio:.10f}',
               '--ratio_max', f'{ratio:.10f}',
               '--n_configs', str(n_configs),
               '--burn_in', str(burn_in),
               '--sweeps_between', str(sweeps_between),
               '--seed', str(seed_base + L*1000 + idx),
               '--ts_file', str(tmp_file)]

        print(f'  [{idx+1:2d}/{len(ratios)}] γ={ratio:.6f} running...', end='', flush=True)
        
        try:
            proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=300)
            if proc.returncode != 0:
                print(f' ERROR ret={proc.returncode}')
                if tmp_file.exists():
                    tmp_file.unlink()
                continue
        except subprocess.TimeoutExpired:
            print(f' TIMEOUT')
            if tmp_file.exists():
                tmp_file.unlink()
            continue
        except Exception as e:
            print(f' FAIL: {e}')
            continue

        # Append rows from tmp_file to out_file
        if tmp_file.exists():
            appended = 0
            try:
                with open(tmp_file) as tf, open(out_file, 'a') as of:
                    for line in tf:
                        if line.strip() == '' or line.startswith('L,'):
                            continue
                        parts = line.split(',')
                        try:
                            if int(parts[0]) != L:
                                continue
                        except Exception:
                            continue
                        of.write(line)
                        appended += 1
                print(f' → {appended} rows')
            except Exception as e:
                print(f' ERROR appending: {e}')
            try:
                tmp_file.unlink()
            except Exception:
                pass
        else:
            print(f' ERROR: tmp file not found')

        time.sleep(0.1)

print('\n'+'='*70)
print('All runs complete. Per-L timeseries files:')
for L in L_vals:
    fpath = bench_dir / f'ts_L{L}_high_res.csv'
    if fpath.exists():
        lines = sum(1 for _ in open(fpath)) - 1  # skip header
        print(f'  {fpath.name}: {lines} rows')
    else:
        print(f'  {fpath.name}: NOT FOUND')
print('='*70)
