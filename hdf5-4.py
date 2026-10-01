
import h5py
import numpy as np
import plotly.graph_objects as go
from datetime import datetime, timezone, timedelta
import os


# ============================================================
# CONFIGURACIÓN
# ============================================================

file_paths = [
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/18_Aug_26/jro20260818_160853.hdf5",
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/19_Aug_26/jro20260819_050853.hdf5",
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/20_Aug_26/jro20260820_050853.hdf5",
    "/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/new-int/clean-test/21_Aug_26/jro20260821_050853.hdf5",
]


# ============================================================
# ARCHIVO DE SALIDA
# ============================================================

directory = os.path.dirname(file_paths[0])

html_file = os.path.join(
    directory,
    "madrigal_interactive_4days.html"
)


# ============================================================
# PARÁMETROS DEL PLOT
# ============================================================

ZMIN = 0
ZMAX = 2.4e12

HEIGHT_MIN = 180
HEIGHT_MAX = 920

N_LEVELS = 15

CONTOUR_SIZE = (
    ZMAX - ZMIN
) / N_LEVELS


# ============================================================
# 1. LEER LOS 4 ARCHIVOS
#
# Primero guardamos los datos de cada archivo.
# Todavía no los concatenamos.
# ============================================================

datasets = []


print()
print("================================================")
print("LEYENDO LOS 4 ARCHIVOS")
print("================================================")


for file_path in file_paths:

    print()
    print("Archivo:")
    print(file_path)

    with h5py.File(file_path, "r") as f:

        ne = f["Data/Array Layout/2D Parameters/ne"][:]

        gdalt = f["Data/Array Layout/gdalt"][:]

        timestamps = f["Data/Array Layout/timestamps"][:]


    print(
        "Shape original de ne:",
        ne.shape
    )

    print(
        "Número de alturas:",
        len(gdalt)
    )

    print(
        "Número de tiempos:",
        len(timestamps)
    )


    # ========================================================
    # 2. ASEGURAR FORMATO
    #
    # ne = altura x tiempo
    # ========================================================

    if ne.shape == (
        len(timestamps),
        len(gdalt)
    ):

        print("Transponiendo ne...")

        ne = ne.T


    elif ne.shape == (
        len(gdalt),
        len(timestamps)
    ):

        print(
            "ne ya está en formato altura x tiempo"
        )


    else:

        raise ValueError(
            f"Dimensiones inesperadas de ne: {ne.shape}"
        )


    print(
        "Shape después de orientación:",
        ne.shape
    )


    # ========================================================
    # 3. SELECCIONAR ALTURAS ENTRE 180 Y 920 km
    # ========================================================

    height_mask = (
        (gdalt >= HEIGHT_MIN) &
        (gdalt <= HEIGHT_MAX)
    )


    gdalt_selected = gdalt[
        height_mask
    ]


    ne = ne[
        height_mask,
        :
    ]


    print(
        "Alturas seleccionadas:",
        gdalt_selected[0],
        "-",
        gdalt_selected[-1],
        "km"
    )


    print(
        "Número de alturas seleccionadas:",
        len(gdalt_selected)
    )


    # ========================================================
    # 4. INTERVALO TEMPORAL NORMAL
    # ========================================================

    dt = np.diff(timestamps)

    dt_normal = np.median(dt)


    print(
        "Intervalo temporal normal:",
        dt_normal,
        "segundos"
    )


    # ========================================================
    # 5. CREAR ÍNDICES TEMPORALES
    # ========================================================

    steps = np.rint(
        (
            timestamps - timestamps[0]
        )
        / dt_normal
    ).astype(int)


    full_steps = np.arange(
        steps.min(),
        steps.max() + 1
    )


    # ========================================================
    # 6. RECONSTRUIR MATRIZ TEMPORAL
    # ========================================================

    ne_full = np.full(
        (
            len(gdalt_selected),
            len(full_steps)
        ),
        np.nan
    )


    indices = (
        steps - steps.min()
    )


    ne_full[
        :,
        indices
    ] = ne


    # ========================================================
    # 7. RECONSTRUIR TIMESTAMPS
    # ========================================================

    timestamps_full = (
        timestamps[0]
        +
        full_steps * dt_normal
    )


    print(
        "Shape final de este día:",
        ne_full.shape
    )

    print(
        "Número de timestamps reconstruidos:",
        len(timestamps_full)
    )


    # ========================================================
    # 8. GUARDAR
    # ========================================================

    datasets.append({

        "ne": ne_full,

        "gdalt": gdalt_selected,

        "timestamps": timestamps_full

    })


# ============================================================
# 9. CREAR LA MALLA VERTICAL COMPLETA
#
# Tomamos TODAS las alturas disponibles en los 4 archivos.
# De esta forma no perdemos las alturas superiores que sólo
# aparecen en algunos archivos.
# ============================================================

all_heights = []


for data in datasets:

    all_heights.extend(
        data["gdalt"]
    )


all_heights = np.array(
    all_heights
)


# ------------------------------------------------------------
# Redondear ligeramente para evitar diferencias numéricas
# como 840.0000001 vs 840.0
# ------------------------------------------------------------

all_heights = np.round(
    all_heights,
    decimals=3
)


# ------------------------------------------------------------
# Obtener alturas únicas
# ------------------------------------------------------------

height_grid = np.unique(
    all_heights
)


height_grid = height_grid[
    (height_grid >= HEIGHT_MIN) &
    (height_grid <= HEIGHT_MAX)
]


height_grid = np.sort(
    height_grid
)


print()
print("================================================")
print("MALLA VERTICAL FINAL")
print("================================================")

print(
    "Altura mínima:",
    height_grid[0],
    "km"
)

print(
    "Altura máxima:",
    height_grid[-1],
    "km"
)

print(
    "Número de alturas:",
    len(height_grid)
)


# ============================================================
# 10. UNIR LOS DATOS TEMPORALES
#
# Primero concatenamos todos los timestamps.
# ============================================================

timestamps_all = np.concatenate(
    [
        data["timestamps"]
        for data in datasets
    ]
)


# ============================================================
# 11. CREAR MATRIZ FINAL
#
# Filas    = alturas de 180 a 920 km
# Columnas = todos los tiempos de los 4 días
#
# Inicialmente TODO es NaN.
# Después colocamos los datos disponibles.
# ============================================================

ne_all = np.full(
    (
        len(height_grid),
        len(timestamps_all)
    ),
    np.nan
)


# ============================================================
# 12. COLOCAR CADA ARCHIVO EN LA MALLA VERTICAL
# ============================================================

column_start = 0


for i, data in enumerate(datasets):

    ne_day = data["ne"]

    gdalt_day = np.round(
        data["gdalt"],
        decimals=3
    )

    n_times = len(
        data["timestamps"]
    )


    column_end = (
        column_start
        +
        n_times
    )


    print()
    print(
        f"Colocando archivo {i + 1}:"
    )

    print(
        "Columnas:",
        column_start,
        "-",
        column_end - 1
    )


    # --------------------------------------------------------
    # Buscar dónde está cada altura del archivo dentro de
    # height_grid
    # --------------------------------------------------------

    for row_day, height in enumerate(
        gdalt_day
    ):

        matches = np.where(
            np.isclose(
                height_grid,
                height,
                atol=0.001
            )
        )[0]


        if len(matches) == 0:

            raise ValueError(
                f"No se pudo ubicar la altura "
                f"{height} km en la malla final."
            )


        row_grid = matches[0]


        ne_all[
            row_grid,
            column_start:column_end
        ] = ne_day[
            row_day,
            :
        ]


    column_start = column_end


# ============================================================
# 13. ORDENAR TODO CRONOLÓGICAMENTE
# ============================================================

sort_indices = np.argsort(
    timestamps_all
)


timestamps_all = (
    timestamps_all[
        sort_indices
    ]
)


ne_all = (
    ne_all[
        :,
        sort_indices
    ]
)


# ============================================================
# 14. ELIMINAR TIMESTAMPS DUPLICADOS
# ============================================================

unique_timestamps, unique_indices = np.unique(
    timestamps_all,
    return_index=True
)


timestamps_all = (
    unique_timestamps
)


ne_all = (
    ne_all[
        :,
        unique_indices
    ]
)


print()
print("================================================")
print("MATRIZ FINAL")
print("================================================")

print(
    "Shape final:",
    ne_all.shape
)

print(
    "Número total de tiempos:",
    len(timestamps_all)
)

print(
    "Número total de alturas:",
    len(height_grid)
)


# ============================================================
# 15. UTC -> LOCAL TIME
#
# Perú = UTC - 5
# ============================================================

times_lt = []


for t in timestamps_all:

    dt_utc = datetime.fromtimestamp(
        t,
        tz=timezone.utc
    )


    dt_local = (
        dt_utc
        -
        timedelta(hours=5)
    )


    times_lt.append(
        dt_local
    )


# ============================================================
# 16. CREAR FIGURA
# ============================================================

fig = go.Figure()


# ============================================================
# 17. CONTOUR
# ============================================================

fig.add_trace(

    go.Contour(

        x=times_lt,

        y=height_grid,

        z=ne_all,


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
# 18. COLOURSCALES DISPONIBLES
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

                    "colorscale": [
                        scale
                    ]

                }

            ]

        )

    )


# ============================================================
# 19. BOTONES SHOW / HIDE LINES
# ============================================================

line_buttons = [

    dict(

        label="Hide lines",

        method="restyle",

        args=[

            {

                "contours.showlines": [
                    False
                ]

            }

        ]

    ),


    dict(

        label="Show lines",

        method="restyle",

        args=[

            {

                "contours.showlines": [
                    True
                ]

            }

        ]

    )

]


# ============================================================
# 20. LAYOUT
# ============================================================

fig.update_layout(

    # --------------------------------------------------------
    # TÍTULO
    # --------------------------------------------------------

    title=dict(

        text="Electron Density — 18–21 August 2026",

        x=0.5

    ),


    # --------------------------------------------------------
    # TAMAÑO
    # --------------------------------------------------------

    width=2000,

    height=700,


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

        range=[
            HEIGHT_MIN,
            HEIGHT_MAX
        ],

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
# 21. TEXTO "COLORSCALE"
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
# 22. TEXTO "LINES"
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
# 23. GUARDAR HTML
# ============================================================

fig.write_html(

    html_file,

    include_plotlyjs=True,

    auto_open=True

)


# ============================================================
# 24. INFORMACIÓN FINAL
# ============================================================

print()
print("==============================================")
print("PLOT INTERACTIVO DE 4 DÍAS GENERADO")
print("==============================================")
print()

print(
    "Archivo:"
)

print(
    html_file
)

print()

print(
    "Altura mínima:",
    HEIGHT_MIN,
    "km"
)

print(
    "Altura máxima del eje:",
    HEIGHT_MAX,
    "km"
)

print(
    "Altura máxima con datos:",
    height_grid[-1],
    "km"
)

print(
    "Número de alturas:",
    len(height_grid)
)

print(
    "Número total de tiempos:",
    len(timestamps_all)
)

print(
    "Shape final:",
    ne_all.shape
)

print()

print(
    "ZMIN:",
    ZMIN
)

print(
    "ZMAX:",
    ZMAX
)

print(
    "Número de niveles:",
    N_LEVELS
)

print(
    "Separación entre niveles:",
    CONTOUR_SIZE
)

print()

print(
    "Colorscale inicial: Jet"
)

print(
    "Lines inicial: Show lines"
)

print()
print("==============================================")

