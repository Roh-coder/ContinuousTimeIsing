#!/usr/bin/env python3
"""
Run benchmark: for each L, run ConTimeIsing to produce a per-config timeseries,
compute sample variance and integrated autocorrelation time of `U`, and
report how many autocorrelation-times are needed to reach error eps=1e-3.

Usage: ./run_benchmark.py [options]

This script requires Python 3 and numpy. matplotlib is optional (for plots).
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from collections import defaultdict

import numpy as np


def parse_ts_file(path):
    # columns: L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site
    rows = []
    with open(path, "r") as f:
        rdr = csv.reader(f)
        try:
            header = next(rdr)
        except StopIteration:
            return rows
        for row in rdr:
            if not row:
                continue
            try:
                L = int(row[0])
                ratio = float(row[1])
                cfg = int(row[2])
                mc_step = int(row[3])
                m2 = float(row[4])
                m4 = float(row[5])
                U = float(row[6])
                avg_segs = float(row[7]) if len(row) > 7 else 0.0
                rows.append({'L': L, 'ratio': ratio, 'cfg': cfg, 'mc_step': mc_step,
                             'm2': m2, 'm4': m4, 'U': U, 'avg_segs': avg_segs})
            except Exception:
                pass
    return rows


def integrated_autocorrelation_time(x, max_lag=None):
    # x: 1D numpy array
    n = len(x)
    if n < 2:
        return 0.0
    x = x - np.mean(x)
    if max_lag is None:
        max_lag = min(n - 1, int(n / 2))

    # FFT-based autocorrelation
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    acf /= acf[0]
    acf = acf[: max_lag + 1]

    # initial positive sequence: sum until first negative autocorr
    tau = 0.5
    for t in range(1, len(acf)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return float(tau)


def run_once(exe, L, n_configs, sweeps_between, burn_in, seed, ts_path):
    cmd = [exe,
           "--layers", str(L),
           "--n_configs", str(n_configs),
           "--sweeps_between", str(sweeps_between),
           "--burn_in", str(burn_in),
           "--seed", str(seed),
           "--ts_file", str(ts_path),
           "--print_every_cfg", "0",
           "--print_every_burn", "0"]

    t0 = time.time()
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Process failed with code {proc.returncode}: {' '.join(cmd)}")
    return elapsed


def analyze(U, eps=1e-3):
    n = len(U)
    if n < 2:
        return None
    mean = float(np.mean(U))
    var = float(np.var(U, ddof=1))
    tau = integrated_autocorrelation_time(U)

    # Required number of configs to reach error eps: N >= (2*tau*var)/eps^2
    N_required = (2.0 * tau * var) / (eps * eps) if eps > 0 else float('inf')
    # In units of autocorrelation times: N_required / (2*tau) = var/eps^2
    ac_units = var / (eps * eps)

    return {
        'n': n,
        'mean': mean,
        'var': var,
        'tau': tau,
        'N_required': N_required,
        'ac_units': ac_units,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--exe', default='./build/ConTimeIsing')
    p.add_argument('--Ls', default='16,32,64,128')
    p.add_argument('--n_configs', type=int, default=500)
    p.add_argument('--sweeps_between', type=int, default=10)
    p.add_argument('--burn_in', type=int, default=200)
    p.add_argument('--seed', type=int, default=1234)
    p.add_argument('--eps', type=float, default=1e-3)
    p.add_argument('--out', default='bench_results.csv')
    args = p.parse_args()

    exe = args.exe
    if not Path(exe).exists():
        print(f"Executable not found: {exe}")
        sys.exit(2)

    Ls = [int(x) for x in args.Ls.split(',') if x.strip()]

    rows = []
    os.makedirs('bench', exist_ok=True)

    for L in Ls:
        print(f"Running L={L} n_configs={args.n_configs} sweeps_between={args.sweeps_between}")
        ts_path = Path('bench') / f'ts_L{L}.csv'
        # run
        elapsed = run_once(exe, L, args.n_configs, args.sweeps_between, args.burn_in, args.seed, ts_path)

        ts_rows = parse_ts_file(ts_path)
        if len(ts_rows) == 0:
            print(f"No timeseries data for L={L} (ts file empty)")
            continue

        # group by ratio
        groups = defaultdict(list)
        for r in ts_rows:
            groups[r['ratio']].append(r['U'])

        total_configs = sum(len(v) for v in groups.values())
        sec_per_config = elapsed / float(total_configs) if total_configs > 0 else float('inf')

        for ratio, Ulist in sorted(groups.items()):
            Uarr = np.array(Ulist, dtype=float)
            info = analyze(Uarr, eps=args.eps)
            if info is None:
                print(f"Not enough data for L={L}, ratio={ratio}")
                continue

            time_required_seconds = info['N_required'] * sec_per_config
            time_required_ac_units = info['N_required'] / (2.0 * info['tau']) if info['tau'] > 0 else float('inf')

            prow = {
                'L': L,
                'ratio': ratio,
                'n_measured': int(info['n']),
                'mean_U': info['mean'],
                'var_U': info['var'],
                'tau_int (configs)': info['tau'],
                'N_required (configs)': info['N_required'],
                'ac_units (N/(2*tau))': info['ac_units'],
                'sec_per_config': sec_per_config,
                'time_required_s': time_required_seconds,
                'time_required_ac_units': time_required_ac_units,
            }
            rows.append(prow)

            print(f"L={L} ratio={ratio} mean_U={info['mean']:.6g} var_U={info['var']:.6g} tau={info['tau']:.3f}")
            print(f" => need {info['N_required']:.1f} configs ({info['ac_units']:.1f} autocorr-times)")
            print(f" elapsed={elapsed:.2f}s sec/config={sec_per_config:.4g} time_required={time_required_seconds:.1f}s")

    keys = ['L','ratio','n_measured','mean_U','var_U','tau_int (configs)','N_required (configs)','ac_units (N/(2*tau))','sec_per_config','time_required_s','time_required_ac_units']
    with open(args.out, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    print('Results written to', args.out)

    # Optionally produce Binder cumulant crossing plot
    try:
        import matplotlib.pyplot as plt
    except Exception:
        print('matplotlib not available; skipping plots')
        return

    # build mapping L -> (ratios, mean_U)
    byL = defaultdict(list)
    for r in rows:
        byL[r['L']].append((r['ratio'], r['mean_U']))

    plt.figure(figsize=(8,6))
    for L, vals in sorted(byL.items()):
        vals_sorted = sorted(vals, key=lambda x: x[0])
        ratios = [v[0] for v in vals_sorted]
        Us = [v[1] for v in vals_sorted]
        plt.plot(ratios, Us, marker='o', label=f'L={L}')

    plt.xlabel('ratio (Gamma)')
    plt.ylabel('Binder U')
    plt.title('Binder cumulant vs ratio')
    plt.legend()
    plt.grid(True)
    plot_path = Path('bench') / 'binder_vs_ratio.png'
    plt.tight_layout()
    plt.savefig(plot_path)
    print('Saved Binder plot to', plot_path)


if __name__ == '__main__':
    main()
