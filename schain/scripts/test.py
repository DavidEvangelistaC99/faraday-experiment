
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

# ==================================
# CONFIGURATION
# ==================================

spectrum_file = "SpectrumDec2024027.dat"
utc_file = "TimeStampsDec2024027.dat"

n_ranges = 300
n_freqs = 200

# ==================================
# READ DATA
# ==================================

# Read spectra
spectra = np.fromfile(spectrum_file, dtype=np.float32)

# Read UTC timestamps
utc = np.fromfile(utc_file, dtype=np.int32)

# Reshape spectra
n_times = len(utc)

spectra = spectra.reshape(n_times, n_ranges, n_freqs)

print("Spectra shape:", spectra.shape)
print("UTC shape:", utc.shape)

# ==================================
# CONVERT UTC TO DATETIME
# ==================================

utc_time = [
    datetime.datetime.fromtimestamp(
        int(t), tz=datetime.timezone.utc
    )
    for t in utc
]

# ==================================
# SELECT RANGE
# ==================================

range_idx = 0

spec = spectra[:, range_idx, :]

# Power in dB
spec_db = 10*np.log10(spec + 1e-12)

# ==================================
# FREQUENCY AXIS
# ==================================

n = 40

freqs = np.fft.fftshift(
    np.fft.fftfreq(n_freqs, d=25e-6)
)#[::2]

print("n_freqs:", n_freqs)
print("freqs:", len(freqs))
print("utc_time:", len(utc_time))
print("spec:", spec.shape)
print("spec_db:", spec_db.shape)
print("spec_db.T:", spec_db.T.shape)

# ==================================
# PLOT
# ==================================

fig, ax = plt.subplots(figsize=(15, 7))

im = ax.pcolormesh(
    utc_time,
    freqs,
    spec_db.T,
    shading="auto",
    cmap="jet",
    vmin=40.0,
    vmax=80.0
)

fig.colorbar(im, ax=ax, label="Power (dB)")

ax.set_xlabel("UTC")
ax.set_ylabel("Frequency (Hz)")

ax.set_title(f"Spectrogram - Range {range_idx}")

ax.xaxis.set_major_formatter(
    mdates.DateFormatter("%H:%M:%S", tz=datetime.timezone.utc)
)

plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

