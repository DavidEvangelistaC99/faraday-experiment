import h5py
import numpy as np
import plotly.graph_objects as go
from datetime import datetime, timezone, timedelta
import os


# ============================================================
# CONFIGURACIÓN
# ============================================================

file_path = (
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/"
    "new-int/clean-test/20_Aug_26/"
    "jro20260820_050853.hdf5"
)

directory = os.path.dirname(file_path)

html_file = os.path.join(
    directory,
    "madrigal_interactive.html"
)


# ============================================================
# PARÁMETROS DEL PLOT
# ============================================================

ZMIN = 0
ZMAX = 2.4e12

# Número aproximado de niveles de contorno
N_LEVELS = 15

# Calcular separación entre niveles
CONTOUR_SIZE = (ZMAX - ZMIN) / N_LEVELS


# ============================================================
# 1. LEER HDF5
# ============================================================

with h5py.File(file_path, "r") as f:

    ne = f["Data/Array Layout/2D Parameters/ne"][:]

    gdalt = f["Data/Array Layout/gdalt"][:]

    timestamps = f["Data/Array Layout/timestamps"][:]


print("Shape original de ne:", ne.shape)
print("Número de alturas:", len(gdalt))
print("Número de tiempos:", len(timestamps))


# ============================================================
# 2. ASEGURAR FORMATO
#
# ne = altura x tiempo
# ============================================================

if ne.shape == (len(timestamps), len(gdalt)):

    print("Transponiendo ne...")

    ne = ne.T

elif ne.shape == (len(gdalt), len(timestamps)):

    print("ne ya está en formato altura x tiempo")

else:

    raise ValueError(
        f"Dimensiones inesperadas de ne: {ne.shape}"
    )


print("Shape final de ne:", ne.shape)


# ============================================================
# 3. INTERVALO TEMPORAL NORMAL
# ============================================================

dt = np.diff(timestamps)

dt_normal = np.median(dt)

print(
    "Intervalo temporal normal:",
    dt_normal,
    "segundos"
)


# ============================================================
# 4. CREAR ÍNDICES TEMPORALES
# ============================================================

steps = np.rint(
    (timestamps - timestamps[0]) / dt_normal
).astype(int)


full_steps = np.arange(
    steps.min(),
    steps.max() + 1
)


# ============================================================
# 5. CREAR MATRIZ COMPLETA
# ============================================================

ne_full = np.full(
    (len(gdalt), len(full_steps)),
    np.nan
)


indices = steps - steps.min()


ne_full[:, indices] = ne


# ============================================================
# 6. RECONSTRUIR TIMESTAMPS
# ============================================================

timestamps_full = (
    timestamps[0]
    + full_steps * dt_normal
)


# ============================================================
# 7. UTC -> LOCAL TIME
#
# Perú = UTC - 5
# ============================================================

times_lt = []

for t in timestamps_full:

    dt_utc = datetime.fromtimestamp(
        t,
        tz=timezone.utc
    )

    dt_local = dt_utc - timedelta(hours=5)

    times_lt.append(dt_local)


# ============================================================
# 8. CREAR FIGURA
# ============================================================

fig = go.Figure()


# ============================================================
# 9. CONTOUR
# ============================================================

fig.add_trace(

    go.Contour(

        x=times_lt,

        y=gdalt,

        z=ne_full,

        # ----------------------------------------------------
        # COLORES
        # ----------------------------------------------------

        colorscale="Jet",

        zmin=ZMIN,

        zmax=ZMAX,

        # ----------------------------------------------------
        # CONTORNOS
        # ----------------------------------------------------

        contours=dict(

            coloring="fill",

            showlines=True,

            start=ZMIN,

            end=ZMAX,

            size=CONTOUR_SIZE
        ),

        # ----------------------------------------------------
        # COLORBAR
        # ----------------------------------------------------

        colorbar=dict(

            title="Ne",

            ticks="outside"
        ),

        # ----------------------------------------------------
        # HOVER
        # ----------------------------------------------------

        hovertemplate=

            "<b>Electron Density</b><br>"
            "Time: %{x|%Y-%m-%d %H:%M:%S}<br>"
            "Height: %{y:.1f} km<br>"
            "Ne: %{z:.3e} m⁻³"
            "<extra></extra>"
    )
)


# ============================================================
# 10. COLOURSCALES DISPONIBLES
# ============================================================

colorscales = [

    "Jet",
    "Viridis",
    "Plasma",
    "Inferno",
    "Turbo",
    "Cividis",
    "Rainbow",
    "Portland",
    "Electric",
    "Earth"
]


colorscale_buttons = []


for scale in colorscales:

    colorscale_buttons.append(

        dict(

            label=scale,

            method="restyle",

            args=[
                {
                    "colorscale": [scale]
                }
            ]
        )
    )


# ============================================================
# 11. BOTONES SHOW / HIDE LINES
# ============================================================

line_buttons = [

    dict(

        label="Hide lines",

        method="restyle",

        args=[
            {
                "contours.showlines": [False]
            }
        ]
    ),

    dict(

        label="Show lines",

        method="restyle",

        args=[
            {
                "contours.showlines": [True]
            }
        ]
    )
]


# ============================================================
# 12. LAYOUT
# ============================================================

fig.update_layout(

    # --------------------------------------------------------
    # TÍTULO
    # --------------------------------------------------------

    title=dict(

        text="Electron Density",

        x=0.5
    ),

    # --------------------------------------------------------
    # TAMAÑO
    # --------------------------------------------------------

    width=1200,

    height=650,

    # --------------------------------------------------------
    # EJE X
    # --------------------------------------------------------

    xaxis=dict(

        title="Time [LT]",

        type="date",

        showgrid=True,

        gridcolor="white"
    ),

    # --------------------------------------------------------
    # EJE Y
    # --------------------------------------------------------

    yaxis=dict(

        title="Range [km]",

        range=[180, 920],

        showgrid=True,

        gridcolor="white"
    ),

    # --------------------------------------------------------
    # HOVER
    # --------------------------------------------------------

    hovermode="closest",

    # --------------------------------------------------------
    # MENÚS
    # --------------------------------------------------------

    updatemenus=[

        # ====================================================
        # COLOURSCALE
        # ====================================================

        dict(

            type="dropdown",

            direction="down",

            x=0.08,

            y=1.08,

            xanchor="left",

            yanchor="top",

            buttons=colorscale_buttons,

            showactive=True,

            bgcolor="white",

            bordercolor="#B0BEC5",

            borderwidth=1
        ),

        # ====================================================
        # LINES
        # ====================================================

        dict(

            type="dropdown",

            direction="down",

            x=0.30,

            y=1.08,

            xanchor="left",

            yanchor="top",

            buttons=line_buttons,

            showactive=True,

            bgcolor="white",

            bordercolor="#B0BEC5",

            borderwidth=1
        )
    ]
)


# ============================================================
# 13. TEXTO "COLORSCALE"
# ============================================================

fig.add_annotation(

    text="Colorscale",

    x=0.02,

    y=1.055,

    xref="paper",

    yref="paper",

    showarrow=False,

    xanchor="left",

    yanchor="middle",

    font=dict(
        size=13
    )
)


# ============================================================
# 14. TEXTO "LINES"
# ============================================================

fig.add_annotation(

    text="Lines",

    x=0.25,

    y=1.055,

    xref="paper",

    yref="paper",

    showarrow=False,

    xanchor="left",

    yanchor="middle",

    font=dict(
        size=13
    )
)


# ============================================================
# 15. GUARDAR HTML
# ============================================================

fig.write_html(

    html_file,

    include_plotlyjs=True,

    auto_open=True
)


# ============================================================
# INFORMACIÓN
# ============================================================

print()
print("==============================================")
print("PLOT INTERACTIVO GENERADO")
print("==============================================")
print()
print("Archivo:")
print(html_file)
print()
print("ZMIN:", ZMIN)
print("ZMAX:", ZMAX)
print("Número de niveles:", N_LEVELS)
print("Separación entre niveles:", CONTOUR_SIZE)
print()
print("Colorscale inicial: Jet")
print("Lines inicial: Show lines")
print()