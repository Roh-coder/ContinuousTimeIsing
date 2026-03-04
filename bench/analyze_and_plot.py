#!/usr/bin/env python3
"""
Analyze timeseries files and generate summary CSV + Binder crossing plot.
"""
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict

def tau_int(x, max_lag=None):
    """Compute integrated autocorrelation time."""
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
    acf = acf[:max_lag + 1]
    
    # Sum until first negative autocorr
    tau = 0.5
    for t in range(1, len(acf)):
        if acf[t] <= 0:
            break
        tau += acf[t]
    return float(tau)

def analyze_ts_files(eps=0.001):
    """Process all ts_L*.csv files and return summary data."""
    results = []
    
    for ts_file in sorted(Path('.').glob('ts_L*.csv')):
        rows = []
        with open(ts_file) as f:
            csv_rdr = csv.reader(f)
            try:
                next(csv_rdr)  # skip header
            except StopIteration:
                continue
            for row in csv_rdr:
                if len(row) > 6:
                    try:
                        L = int(row[0])
                        ratio = float(row[1])
                        U = float(row[6])
                        rows.append({'L': L, 'ratio': ratio, 'U': U})
                    except (ValueError, IndexError):
                        pass
        
        if not rows:
            continue
        
        # Group by ratio
        groups = defaultdict(list)
        for r in rows:
            groups[r['ratio']].append(r['U'])
        
        L = rows[0]['L']
        for ratio in sorted(groups.keys()):
            Uarr = np.array(groups[ratio])
            mean_U = float(np.mean(Uarr))
            var_U = float(np.var(Uarr, ddof=1))
            tau = tau_int(Uarr)
            N_required = (2.0 * tau * var_U) / (eps * eps) if eps > 0 else float('inf')
            ac_units = var_U / (eps * eps)
            
            results.append({
                'L': L,
                'ratio': ratio,
                'n_measured': len(Uarr),
                'mean_U': mean_U,
                'var_U': var_U,
                'tau_int (configs)': tau,
                'N_required (configs)': N_required,
                'ac_units (N/(2*tau))': ac_units,
            })
    
    return results

def write_summary_csv(results, filename='summary_per_point.csv'):
    """Write summary to CSV."""
    if not results:
        print("No results to write")
        return
    
    keys = ['L', 'ratio', 'n_measured', 'mean_U', 'var_U', 
            'tau_int (configs)', 'N_required (configs)', 'ac_units (N/(2*tau))']
    
    with open(filename, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in sorted(results, key=lambda x: (x['L'], x['ratio'])):
            wr.writerow(r)
    
    print(f"Summary written to {filename}")
    return filename

def make_binder_plot(results, filename='binder_vs_ratio.png'):
    """Create Binder U vs ratio crossing plot."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return None
    
    # Group by L
    byL = defaultdict(list)
    for r in results:
        byL[r['L']].append((r['ratio'], r['mean_U']))
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for L in sorted(byL.keys()):
        vals = sorted(byL[L], key=lambda x: x[0])
        ratios = [v[0] for v in vals]
        Us = [v[1] for v in vals]
        ax.plot(ratios, Us, marker='o', linewidth=2, markersize=8, label=f'L={L}')
    
    ax.set_xlabel('Coupling ratio (Gamma)', fontsize=12)
    ax.set_ylabel('Binder Cumulant U', fontsize=12)
    ax.set_title('Binder Cumulant Crossing Diagram', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Plot saved to {filename}")
    return filename

def main():
    print("Analyzing timeseries files...")
    results = analyze_ts_files(eps=0.001)
    
    if not results:
        print("No data found")
        return
    
    print(f"Found {len(results)} data points")
    
    # Print summary to screen
    print("\n" + "="*90)
    print(f"{'L':>4} {'ratio':>8} {'n':>6} {'mean_U':>10} {'var_U':>12} {'tau':>10} {'ac_units':>12}")
    print("="*90)
    for r in sorted(results, key=lambda x: (x['L'], x['ratio'])):
        print(f"{r['L']:>4d} {r['ratio']:>8.2f} {r['n_measured']:>6d} {r['mean_U']:>10.6f} "
              f"{r['var_U']:>12.6f} {r['tau_int (configs)']:>10.2f} {r['ac_units (N/(2*tau))']:>12.1f}")
    print("="*90)
    
    # Write CSV
    write_summary_csv(results)
    
    # Make plot
    make_binder_plot(results)

if __name__ == '__main__':
    main()
