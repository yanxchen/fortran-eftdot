import numpy as np
import accupy
import struct
import argparse

def generate_binary_data(n=100, samples=1000):
    target_conds = 10**np.linspace(0, 35, 20)
    samples_per_tc = samples // len(target_conds)

    with open('input_data.bin', 'wb') as f:
        for tc in target_conds:
            for _ in range(samples_per_tc):
                # Generate ill-conditioned vectors
                x, y, exact, actual_cond = accupy.generate_ill_conditioned_dot_product(n, tc)
                
                # Pack metadata (as DP)
                f.write(struct.pack('d', float(actual_cond)))
                f.write(struct.pack('d', float(exact)))
                
                # Pack SP vectors (binary32) - Directly read into real(sp)
                f.write(struct.pack(f'{n}f', *x.astype(np.float32)))
                f.write(struct.pack(f'{n}f', *y.astype(np.float32)))
                
                # Pack DP vectors (binary64) - Directly read into real(dp)
                f.write(struct.pack(f'{n}d', *x.astype(np.float64)))
                f.write(struct.pack(f'{n}d', *y.astype(np.float64)))
                
    print(f"Generated {samples} samples with dedicated SP and DP vector sets.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate Binary Test Data for EFTDot')
    parser.add_argument('--n', type=int, default=100, help='Vector length')
    parser.add_argument('--samples', type=int, default=1000, help='Number of samples')
    args = parser.parse_args()

    n = args.n
    samples = args.samples

    generate_binary_data(n, samples)