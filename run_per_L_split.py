#!/usr/bin/env python3
"""
Run ConTimeIsing across many ratios, saving each L into its own timeseries CSV.
- Resumes if a per-L CSV already contains a ratio.
- Writes tmp per-job file and appends to per-L file on success.
Usage: python3 run_per_L_split.py
"""
import subprocess
import sys
import time
from pathlib import Path
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--L', nargs='+', type=int,
                    default=[32,28,24,20,16,12,8,4],
                    help='List of L values (will run in this order)')
parser.add_argument('--ratio_min', type=float, default=0.9)
parser.add_argument('--ratio_max', type=float, default=10.0/9.0)
parser.add_argument('--n_ratios', type=int, default=21)
parser.add_argument('--n_configs', type=int, default=500)
parser.add_argument('--burn_in', type=int, default=250)
parser.add_argument('--sweeps_between', type=int, default=5)
parser.add_argument('--seed_base', type=int, default=54321)
parser.add_argument('--concurrency_delay', type=float, default=0.1,
                    help='Delay between jobs (s)')
parser.add_argument('--dry_run', action='store_true')
args = parser.parse_args()

L_vals = args.L
ratios = list(reversed([args.ratio_min + i*(args.ratio_max-args.ratio_min)/(args.n_ratios-1)
                        for i in range(args.n_ratios)]))

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
                    existing_ratios.add(float(r['ratio']))
                except Exception:
                    pass
    else:
        # create header
        with open(out_file, 'w') as f:
            f.write('L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n')

    print(f'Running L={L} → existing ratios: {sorted(existing_ratios)}')

    for idx, ratio in enumerate(ratios):
        if ratio in existing_ratios:
            print(f'  Skipping ratio={ratio:.6f} (already present)')
            continue

        tmp_file = bench_dir / f'tmp_ts_L{L}_r{ratio:.6f}.csv'
        cmd = [str(con_exe),
               '--layers', str(L),
               '--ratio_min', f'{ratio:.10f}',
               '--ratio_max', f'{ratio:.10f}',
               '--n_configs', str(args.n_configs),
               '--burn_in', str(args.burn_in),
               '--sweeps_between', str(args.sweeps_between),
               '--seed', str(args.seed_base + L*1000 + idx),
               '--ts_file', str(tmp_file)]

        print(f'  Running ratio={ratio:.6f} (seed={args.seed_base + L*1000 + idx})')
        if args.dry_run:
            print('   dry-run: would execute:', ' '.join(cmd))
            continue

        # run the simulator
        try:
            proc = subprocess.run(cmd, check=False)
            if proc.returncode != 0:
                print(f'   ERROR: process returned {proc.returncode}; skipping this ratio')
                if tmp_file.exists():
                    tmp_file.unlink()
                continue
        except FileNotFoundError:
            print('   ERROR: ConTimeIsing executable not found or not executable')
            sys.exit(1)

        # Append tmp_file rows for this L to out_file (skip header)
        if tmp_file.exists():
            appended = 0
            with open(tmp_file) as tf, open(out_file, 'a') as of:
                for line in tf:
                    if line.strip() == '':
                        continue
                    if line.startswith('L,'):
                        continue
                    # only append rows matching this L (defensive)
                    parts = line.split(',')
                    try:
                        if int(parts[0]) != L:
                            continue
                    except Exception:
                        continue
                    of.write(line)
                    appended += 1
            print(f'   Appended {appended} rows to {out_file}')
            try:
                tmp_file.unlink()
            except Exception:
                pass
        else:
            print('   WARNING: expected tmp file not found; nothing appended')

        time.sleep(args.concurrency_delay)

print('\nAll done. Per-L time series files are in bench/ (ts_L<...>_high_res.csv)')
