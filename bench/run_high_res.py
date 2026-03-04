#!/usr/bin/env python3
"""
High-resolution benchmark around critical point.
- Tight ratio range (0.85-1.15, 13 points)
- Multiple lattice sizes (16, 20, 24, 28, 32, 40, 48, 56, 64)
- High accuracy (500 configs per point)
"""
import subprocess
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict


def parse_ts_file(path):
    """Parse timeseries CSV: L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site"""
    rows = []
    with open(path, "r") as f:
        rdr = csv.reader(f)
        try:
            next(rdr)  # skip header
        except StopIteration:
            return rows
        for row in rdr:
            if not row:
                continue
            try:
                rows.append({
                    'L': int(row[0]),
                    'ratio': float(row[1]),
                    'U': float(row[6])
                })
            except Exception:
                pass
    return rows


def integrated_autocorrelation_time(x):
    """FFT-based autocorrelation time."""
    n = len(x)
    if n < 2:
        return 0.5
    x = x - np.mean(x)
    f = np.fft.rfft(np.concatenate([x, np.zeros_like(x)]))
    acf = np.fft.irfft(f * np.conjugate(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for t in range(1, min(len(acf), n // 2)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return tau


def run_benchmark_high_res():
    """Run high-resolution benchmark around critical point."""
    
    # Tight ratio range around criticality
    ratios = np.linspace(0.85, 1.15, 13)
    ratio_str = ",".join([f"{r:.4f}" for r in ratios])
    
    # Multiple lattice sizes for better finite-size scaling
    L_values = [16, 20, 24, 28, 32, 40, 48, 56, 64]
    L_str = ",".join([str(L) for L in L_values])
    
    # High accuracy: 500 configs per point
    n_configs = 500
    sweeps_between = 5
    burn_in = 200
    
    # Output file
    ts_file = "bench/ts_high_res.csv"
    summary_file = "bench/summary_high_res.csv"
    
    print(f"Running high-resolution benchmark:")
    print(f"  Ratios: {len(ratios)} points from {ratios[0]:.3f} to {ratios[-1]:.3f}")
    print(f"  L values: {L_values}")
    print(f"  Configs per point: {n_configs}")
    print(f"  Output: {ts_file}")
    print()
    
    # Build command
    cmd = [
        "./build/ConTimeIsing",
        "--layers", L_str,
        "--ratio_min", str(ratios[0]),
        "--ratio_max", str(ratios[-1]),
        "--ratio_n", str(len(ratios)),
        "--n_configs", str(n_configs),
        "--sweeps_between", str(sweeps_between),
        "--burn_in", str(burn_in),
        "--seed", "42",
        "--BC", "1",
        "--ts_file", ts_file,
        "--print_every_burn", "100",
        "--print_every_cfg", "50"
    ]
    
    print("Running: " + " ".join(cmd))
    print()
    
    # Run simulator
    result = subprocess.run(cmd, cwd="/workspaces/ContinuousTimeIsing", capture_output=False)
    if result.returncode != 0:
        print(f"Error: simulator exit code {result.returncode}")
        return False
    
    # Parse results
    print("\nAnalyzing results...")
    rows = parse_ts_file(ts_file)
    print(f"Read {len(rows)} data points from {ts_file}")
    
    # Group by (L, ratio)
    eps = 0.001
    results = []
    
    by_L_ratio = defaultdict(list)
    for r in rows:
        key = (r['L'], round(r['ratio'], 6))
        by_L_ratio[key].append(r['U'])
    
    for (L, ratio), U_list in sorted(by_L_ratio.items()):
        U_arr = np.array(U_list)
        mean_U = np.mean(U_arr)
        var_U = np.var(U_arr, ddof=1)
        n_measured = len(U_arr)
        
        # Autocorrelation time
        tau = integrated_autocorrelation_time(U_arr)
        
        # Cost metrics
        ac_units = var_U / (eps ** 2)
        N_req = 2 * tau * var_U / (eps ** 2)
        
        results.append({
            'L': L,
            'ratio': ratio,
            'n_measured': n_measured,
            'mean_U': mean_U,
            'var_U': var_U,
            'tau_int': tau,
            'ac_units': ac_units,
            'N_required': N_req
        })
        
        print(f"  L={L:2d}, ratio={ratio:.4f}: n={n_measured}, U={mean_U:.5f}, var={var_U:.6f}, tau={tau:.2f}, ac_u={ac_units:.0f}")
    
    # Write summary
    with open(summary_file, 'w') as f:
        f.write('L,ratio,n_measured,mean_U,var_U,tau_int (configs),N_required (configs),ac_units (N/(2*tau))\n')
        for r in sorted(results, key=lambda x: (x['L'], x['ratio'])):
            f.write(f"{r['L']},{r['ratio']},{r['n_measured']},{r['mean_U']:.17g},"
                   f"{r['var_U']:.17g},{r['tau_int']:.6f},{r['N_required']:.2f},"
                   f"{r['ac_units']:.1f}\n")
    
    print(f"\nWrote summary to {summary_file}")
    return True


if __name__ == "__main__":
    run_benchmark_high_res()
