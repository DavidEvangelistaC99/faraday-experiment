import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from datetime import datetime, timezone, timedelta


# =========================================================
# 1. ARCHIVOS
# =========================================================

files = [
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/18_Aug_26/jro20260818_160853.hdf5",

    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/19_Aug_26/jro20260819_050853.hdf5",

    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/20_Aug_26/jro20260820_050853.hdf5",

    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/21_Aug_26/jro20260821_050853.hdf5"
]


# =========================================================
# 2. LEER LOS ARCHIVOS
# =========================================================

data = []

all_heights = []


for filename in files:

    print("\nLeyendo:")
    print(filename)

    with h5py.File(filename, "r") as f:

        ne = f["Data/Array Layout/2D Parameters/ne"][:]

        gdalt = f["Data/Array Layout/gdalt"][:]

        timestamps = f["Data/Array Layout/timestamps"][:]


    print("ne:", ne.shape)
    print("gdalt:", gdalt.shape)
    print("timestamps:", timestamps.shape)


    # -----------------------------------------------------
    # Verificar orientación
    # -----------------------------------------------------

    if ne.shape == (len(timestamps), len(gdalt)):
        ne = ne.T

    elif ne.shape != (len(gdalt), len(timestamps)):

        raise ValueError(
            "Las dimensiones de ne no corresponden con "
            "gdalt y timestamps."
        )


    data.append({
        "ne": ne,
        "gdalt": gdalt,
        "timestamps": timestamps
    })


    #all_heights.extend(gdalt)
    all_heights.append(gdalt)


# =========================================================
# 3. CREAR VECTOR DE ALTURAS COMÚN
# =========================================================

gdalt_common = np.unique(
    np.concatenate(all_heights)
)


print("\nAlturas comunes:", len(gdalt_common))


# =========================================================
# 4. CREAR MATRICES CON ALTURAS COMUNES
# =========================================================

ne_list = []
timestamps_list = []


for d in data:

    ne = d["ne"]
    gdalt = d["gdalt"]
    timestamps = d["timestamps"]


    # Matriz llena de NaN

    ne_common = np.full(
        (len(gdalt_common), len(timestamps)),
        np.nan
    )


    # Buscar dónde está cada altura del archivo
    # dentro del vector común

    for i, h in enumerate(gdalt):

        idx = np.argmin(
            np.abs(gdalt_common - h)
        )

        ne_common[idx, :] = ne[i, :]


    ne_list.append(ne_common)

    timestamps_list.append(timestamps)


# =========================================================
# 5. CONCATENAR LOS CUATRO ARCHIVOS
# =========================================================

ne_all = np.concatenate(
    ne_list,
    axis=1
)


timestamps_all = np.concatenate(
    timestamps_list
)


print("\nDatos unidos:")
print("ne_all:", ne_all.shape)
print("timestamps_all:", timestamps_all.shape)
print("gdalt_common:", gdalt_common.shape)


# =========================================================
# 6. ORDENAR CRONOLÓGICAMENTE
# =========================================================

order = np.argsort(timestamps_all)

timestamps_all = timestamps_all[order]

ne_all = ne_all[:, order]


# =========================================================
# 7. ELIMINAR TIMESTAMPS DUPLICADOS
# =========================================================

timestamps_unique, unique_indices = np.unique(
    timestamps_all,
    return_index=True
)

ne_all = ne_all[:, unique_indices]

timestamps_all = timestamps_unique


# =========================================================
# 8. CALCULAR INTERVALO NORMAL
# =========================================================

dt = np.diff(timestamps_all)

dt_normal = np.median(dt)

print("\nIntervalo normal:")
print(dt_normal, "segundos")
print(dt_normal / 60, "minutos")


# =========================================================
# 9. CREAR GRILLA TEMPORAL COMPLETA
# =========================================================

steps = np.rint(
    (timestamps_all - timestamps_all[0]) / dt_normal
).astype(int)


full_steps = np.arange(
    steps[0],
    steps[-1] + 1
)


# =========================================================
# 10. MATRIZ FINAL CON NaN
# =========================================================

ne_full = np.full(
    (len(gdalt_common), len(full_steps)),
    np.nan
)


indices = steps - steps[0]

ne_full[:, indices] = ne_all


# =========================================================
# 11. TIMESTAMPS COMPLETOS
# =========================================================

timestamps_full = (
    timestamps_all[0]
    + full_steps * dt_normal
)


# =========================================================
# 12. UTC -> LOCAL TIME
# =========================================================

times_lt = np.array([

    datetime.fromtimestamp(
        t,
        tz=timezone.utc
    ) - timedelta(hours=5)

    for t in timestamps_full

])


# =========================================================
# 13. MATPLOTLIB
# =========================================================

times_num = mdates.date2num(times_lt)


# =========================================================
# 14. PLOT
# =========================================================

fig, ax = plt.subplots(
    figsize=(16, 3.5)
)


'''pcm = ax.pcolormesh(

    times_num,

    gdalt_common,

    ne_full,

    shading="nearest",

    cmap="jet"

)'''


pcm = ax.pcolormesh(
    times_num,
    gdalt_common,
    ne_full,
    cmap='jet',
    vmin=0,
    vmax=2e12,#,
    shading='nearest'
)


levels = np.linspace(
    np.nanmin(ne_full),
    np.nanmax(ne_full),
    10
)

cs = ax.contour(
    times_num,
    gdalt_common,
    ne_full,
    levels=levels,
    colors='black',
    linewidths=0.6,
    alpha=0.5
)
# =========================================================
# 15. EJES
# =========================================================

ax.set_xlabel(
    "Time (LT)"
)

ax.set_ylabel(
    "Range (km)"
)
ax.set_ylim(180, 920)


ax.xaxis_date()

ax.xaxis.set_major_formatter(
    mdates.DateFormatter("%d-%m %H:%M")
)


# =========================================================
# 16. COLORBAR
# =========================================================

cbar = fig.colorbar(
    pcm,
    ax=ax
)

cbar.set_label(
    "Ne"
)


# =========================================================
# 17. MOSTRAR
# =========================================================

plt.tight_layout()
plt.savefig('Ne.png')
plt.show()