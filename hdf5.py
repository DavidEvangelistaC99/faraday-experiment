import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timezone, timedelta
from matplotlib.colors import LogNorm
import os

#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/clean-test/18_Aug_26-2/jro20260818_161200.hdf5"
#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/clean-test/19_Aug_26-5/jro20260819_050001.hdf5"
#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/clean-test/20_Aug_26-6/jro20260820_050001.hdf5"
#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/clean-test/21_Aug_26-3/jro20260821_051201.hdf5"

#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/18_Aug_26/jro20260818_160853.hdf5"

#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/19_Aug_26/jro20260819_050853.hdf5"
#filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/20_Aug_26/jro20260820_050853.hdf5"
# filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/21_Aug_26/jro20260821_050853.hdf5"

filename = "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/20_Aug_26/jro20260820_050853.hdf5"

directory = os.path.dirname(filename)
# ---------------------------------------------------------
# 1. Leer HDF5
# ---------------------------------------------------------

with h5py.File(filename, "r") as f:

    ne = f["Data/Array Layout/2D Parameters/ne"][:]
    gdalt = f["Data/Array Layout/gdalt"][:]
    timestamps = f["Data/Array Layout/timestamps"][:]


# ---------------------------------------------------------
# 2. Verificar dimensiones
# ---------------------------------------------------------

print("ne:", ne.shape)
print("gdalt:", gdalt.shape)
print("timestamps:", timestamps.shape)


# ne debe ser:
#
#       altura x tiempo
#
#       (48, 58)


if ne.shape == (len(timestamps), len(gdalt)):
    ne = ne.T


# ---------------------------------------------------------
# 3. Encontrar el intervalo temporal normal
# ---------------------------------------------------------

dt = np.diff(timestamps)

dt_normal = np.median(dt)

print("Intervalo normal:", dt_normal, "s")


# ---------------------------------------------------------
# 4. Convertir cada timestamp a un índice temporal
# ---------------------------------------------------------

steps = np.rint(
    (timestamps - timestamps[0]) / dt_normal
).astype(int)


# ---------------------------------------------------------
# 5. Crear grilla temporal completa
# ---------------------------------------------------------

full_steps = np.arange(
    steps[0],
    steps[-1] + 1
)


# Matriz completa inicialmente llena de NaN

ne_full = np.full(
    (len(gdalt), len(full_steps)),
    np.nan
)


# ---------------------------------------------------------
# 6. Colocar los datos originales
# ---------------------------------------------------------

indices = steps - steps[0]

ne_full[:, indices] = ne


# ---------------------------------------------------------
# 7. Crear timestamps completos
# ---------------------------------------------------------

timestamps_full = (
    timestamps[0]
    + full_steps * dt_normal
)


# ---------------------------------------------------------
# 8. Convertir UTC -> LT
# ---------------------------------------------------------

times_lt = np.array([
    datetime.fromtimestamp(
        t,
        tz=timezone.utc
    ) - timedelta(hours=5)
    for t in timestamps_full
])


# ---------------------------------------------------------
# 9. Convertir fechas a formato matplotlib
# ---------------------------------------------------------

times_num = mdates.date2num(times_lt)


# ---------------------------------------------------------
# 10. PLOT
# ---------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(12, 3.5)
)


'''pcm = ax.pcolormesh(
    times_num,
    gdalt,
    ne_full,
    cmap='jet'
)'''

pcm = ax.pcolormesh(
    times_num,
    gdalt,
    ne_full,
    cmap='jet',
    shading='nearest',
    vmin=0,
    vmax=2e12
)

'''pcm = ax.pcolormesh(
    times_num,
    gdalt,
    ne_full,
    shading="nearest",
    cmap="jet",
    norm=LogNorm(vmin=1e4, vmax=1e7)
)'''


ax.set_xlabel("Time (LT)")
ax.set_ylabel("Range (km)")
ax.set_ylim(180, 920)

ax.xaxis_date()

ax.xaxis.set_major_formatter(
    mdates.DateFormatter("%H:%M")
)


fig.colorbar(
    pcm,
    ax=ax,
    label="Ne"
)


plt.tight_layout()
plt.savefig(directory+'/madrigal.png')
plt.show()