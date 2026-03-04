#!/usr/bin/env python3
"""
Run targeted benchmark: L=48,32,24,16 with 15 ratios (10/9 to 0.9), 500 configs each.
Resumes by skipping (L,ratio) pairs already present in per-L CSV files.
"""
import subprocess
import sys
import time
from pathlib import Path
import csv
import numpy as np

L_vals = [48, 32, 24, 16]
ratios = np.linspace(10.0/9.0, 0.9, 15)
n_configs = 500
burn_in = 250
sweeps_between = 5
seed_base = 99999

bench_dir = Path('bench')
bench_dir.mkdir(exist_ok=True)

con_exe = Path('build/ConTimeIsing')
if not con_exe.exists():
    print('Error: build/ConTimeIsing not found')
    sys.exit(1)

total_jobs = len(L_vals) * len(ratios)
completed = 0

for L in L_vals:
    out_file = bench_dir / f'ts_L{L}_high_res.csv'
    existing_ratios = set()
    
    # Read existing ratios in this L's file
    if out_file.exists():
        with open(out_file) as f:
            reader = csv.DictReader(f)
            for r in reader:
                try:
                    existing_ratios.add(round(float(r['ratio']), 6))
                except Exception:
                    pass
    else:
        # Create header
        with open(out_file, 'w') as f:
            f.write('L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n')

    print(f'\nL={L}: {len(existing_ratios)} ratios already present')

    for idx, ratio in enumerate(ratios):
        r_rounded = round(ratio, 6)
        if r_rounded in existing_ratios:
            print(f'  Skipping γ={ratio:.6f} (present)')
            completed += 1
            continue

        tmp_file = bench_dir / f'tmp_L{L}_r{ratio:.8f}.csv'
        cmd = [str(con_exe),
               '--layers', str(L),
               '--ratio_min', f'{ratio:.10f}',
               '--ratio_max', f'{ratio:.10f}',
               '--n_configs', str(n_configs),
               '--burn_in', str(burn_in),
               '--sweeps_between', str(sweeps_between),
               '--seed', str(seed_base + L*1000 + idx),
               '--ts_file', str(tmp_file)]

        print(f'  Running γ={ratio:.6f} [{completed+1}/{total_jobs}]')
        try:
            proc = subprocess.run(cmd, check=False, timeout=600)
            if proc.returncode != 0:
                print(f'    ERROR: returned {proc.returncode}')
                if tmp_file.exists():
                    tmp_file.unlink()
                continue
        except subprocess.TimeoutExpired:
            print(f'    TIMEOUT after 600s')
            if tmp_file.exists():
                tmp_file.unlink()
            continue
        except Exception as e:
            print(f'    ERROR: {e}')
            if tmp_file.exists():
                tmp_file.unlink()
            continue

        # Append to per-L file
        if tmp_file.exists():
            rows_added = 0
            try:
                with open(tmp_file) as tf, open(out_file, 'a') as of:
                    for line in tf:
                        if line.strip() == '' or line.startswith('L,'):
                            continue
                        of.write(line)
                        rows_added += 1
                print(f'    Appended {rows_added} rows')
            except Exception as e:
                print(f'    ERROR appending: {e}')
            finally:
                try:
                    tmp_file.unlink()
                except Exception:
                    pass
        
        completed += 1
        time.sleep(0.2)

print(f'\nDone: {completed}/{total_jobs} jobs completed')
