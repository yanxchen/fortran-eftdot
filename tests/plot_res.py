import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_comprehensive_benchmark(csv_file='results.csv'):
    # 1. Load and clean the data
    # 'index_col=False' and 'dropna' handle the potential trailing comma issue
    df = pd.read_csv(csv_file, index_col=False).dropna(axis=1, how='all')
    
    # Convert all columns to numeric, ensuring the 10^36 values are read correctly
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna().query('cond > 0')

    print(f"Plotting {len(df)} samples.")
    print(f"Condition Range: {df['cond'].min():.1e} to {df['cond'].max():.1e}")

    # 2. Setup constants
    n = 100
    eps_s, eps_d = 2**-24, 2**-53
    
    def get_rel_err(calculated, reference):
        with np.errstate(divide='ignore', invalid='ignore'):
            err = np.abs(calculated - reference) / np.abs(reference)
        # Cap at 2.0 to show the failure plateau clearly
        return np.where(np.isfinite(err), np.minimum(err, 2.0), 2.0)

    plt.figure(figsize=(12, 8))

    # 3. Plot Single Precision (SP) Group
    # Dot2 (SP) - Blue dots
    plt.scatter(df['cond'], get_rel_err(df['dot3_s'], df['accurate']), 
                s=15, color='blue', alpha=0.6, label='Dot2 (SP)', zorder=5)
    
    # Intrinsic (SP) - Cyan dots
    plt.scatter(df['cond'], get_rel_err(df['int_s'], df['accurate']), 
                s=8, color='cyan', alpha=0.3, label='Intrinsic (SP)', zorder=3)
    
    # Manual Loop (SP) - Small 'x' markers
    plt.scatter(df['cond'], get_rel_err(df['loop_s'], df['accurate']), 
                s=10, marker='x', color='black', alpha=0.3, label='Loop (SP)', zorder=2)

    # 4. Plot Double Precision (DP) Group
    # Dot2 (DP) - Red dots
    plt.scatter(df['cond'], get_rel_err(df['dot3_d'], df['accurate']), 
                s=15, color='red', alpha=0.6, label='Dot2 (DP)', zorder=5)
    
    # Intrinsic (DP) - Orange dots
    plt.scatter(df['cond'], get_rel_err(df['int_d'], df['accurate']), 
                s=8, color='orange', alpha=0.3, label='Intrinsic (DP)', zorder=3)
    
    # Manual Loop (DP) - Small '+' markers
    plt.scatter(df['cond'], get_rel_err(df['loop_d'], df['accurate']), 
                s=10, marker='+', color='brown', alpha=0.3, label='Loop (DP)', zorder=2)

    # 5. Plot Theoretical Bounds (Formula 5.8)
    c_axis = np.logspace(0, 36, 1000)
    g_s, g_d = (n*eps_s)/(1-n*eps_s), (n*eps_d)/(1-n*eps_d)
    
    plt.plot(c_axis, np.minimum(eps_s + 0.5 * g_s**2 * c_axis, 2.0), 
             'b--', linewidth=1.5, label='Theory (SP)')
    plt.plot(c_axis, np.minimum(eps_d + 0.5 * g_d**2 * c_axis, 2.0), 
             'r--', linewidth=1.5, label='Theory (DP)')

    # 6. Formatting (CRITICAL for log-log visualization)
    plt.xscale('log')
    plt.yscale('log')
    
    # Force the plot to show the full span of 36 orders of magnitude
    plt.xlim(1e0, 1e36)
    plt.ylim(1e-19, 10)

    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.xlabel('Condition Number $cond(x^T y)$')
    plt.ylabel('Relative Error')
    plt.title('Dot products Errors')
    plt.legend(loc='lower right', ncol=2, fontsize='small')
    
    plt.tight_layout()
    plt.savefig('full_comparison_plot.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_comprehensive_benchmark()