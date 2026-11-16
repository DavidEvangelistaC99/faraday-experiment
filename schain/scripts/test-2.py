import numpy as np

spectra = np.fromfile("SpectrumDec2024027.dat", dtype=np.float32)
utc = np.fromfile("TimeStampsDec2024027.dat", dtype=np.int32)

n_ranges = 300
n_freqs = 40

n_spectra = len(spectra) // (n_ranges*n_freqs)

print("Cantidad de UTC:", len(utc))
print("Cantidad de espectros:", n_spectra)
print("Elementos sobrantes:", len(spectra) % (n_ranges*n_freqs))