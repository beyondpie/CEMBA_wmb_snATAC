from typing import List
import math

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# * matplot color introduction
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# https://matplotlib.org/stable/tutorials/colors/colormapnorms.html

# * color palettable
# https://jiffyclub.github.io/palettable/
# https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html

def plot_colortable(colors, *, ncols=4, sort_colors=True):
    """Show colors.

    Ref: https://matplotlib.org/stable/gallery/color/named_colors.html
    Usage:
    plot_colortable(mcolors.CSS4_COLORS)
    plt.show()
    """

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12

    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        names = sorted(
            colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c)))
        )
    else:
        names = list(colors)

    n = len(names)
    nrows = math.ceil(n / ncols)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(
        margin / width,
        margin / height,
        (width - margin) / width,
        (height - margin) / height,
    )
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows - 0.5), -cell_height / 2.0)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(
            text_pos_x,
            y,
            name,
            fontsize=14,
            horizontalalignment="left",
            verticalalignment="center",
        )

        ax.add_patch(
            Rectangle(
                xy=(swatch_start_x, y - 9),
                width=swatch_width,
                height=18,
                facecolor=colors[name],
                edgecolor="0.7",
            )
        )

    return fig


# colors from SnapATAC
SnapATACPalette: List[str] = [
    "#E31A1C",
    "#FFD700",
    "#771122",
    "#777711",
    "#1F78B4",
    "#68228B",
    "#AAAA44",
    "#60CC52",
    "#771155",
    "#DDDD77",
    "#774411",
    "#AA7744",
    "#AA4455",
    "#117744",
    "#000080",
    "#44AA77",
    "#AA4488",
    "#DDAA77",
    "#D9D9D9",
    "#BC80BD",
    "#FFED6F",
    "#7FC97F",
    "#BEAED4",
    "#FDC086",
    "#FFFF99",
    "#386CB0",
    "#F0027F",
    "#BF5B17",
    "#666666",
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#A6CEE3",
    "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FB9A99",
    "#E31A1C",
    "#FDBF6F",
    "#FF7F00",
    "#CAB2D6",
    "#6A3D9A",
    "#B15928",
    "#FBB4AE",
    "#B3CDE3",
    "#CCEBC5",
    "#DECBE4",
    "#FED9A6",
    "#FFFFCC",
    "#E5D8BD",
    "#FDDAEC",
    "#F2F2F2",
    "#B3E2CD",
    "#FDCDAC",
    "#CBD5E8",
    "#F4CAE4",
    "#E6F5C9",
    "#FFF2AE",
    "#F1E2CC",
    "#CCCCCC",
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FFFF33",
    "#A65628",
    "#F781BF",
    "#999999",
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3",
    "#8DD3C7",
    "#FFFFB3",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#FCCDE5",
]

# colors from ArchR
stallion: List[str] = [
    "#D51F26",
    "#272E6A",
    "#208A42",
    "#89288F",
    "#F47D2B",
    "#FEE500",
    "#8A9FD1",
    "#C06CAB",
    "#E6C2DC",
    "#90D5E4",
    "#89C75F",
    "#F37B7D",
    "#9983BD",
    "#D24B27",
    "#3BBCA8",
    "#6E4B9E",
    "#0C727C",
    "#7E1416",
    "#D8A767",
    "#3D3D3D",
]

stallion2: List[str] = [
    "#D51F26",
    "#272E6A",
    "#208A42",
    "#89288F",
    "#F47D2B",
    "#FEE500",
    "#8A9FD1",
    "#C06CAB",
    "#E6C2DC",
    "#90D5E4",
    "#89C75F",
    "#F37B7D",
    "#9983BD",
    "#D24B27",
    "#3BBCA8",
    "#6E4B9E",
    "#0C727C",
    "#7E1416",
    "#D8A767",
]

calm: List[str] = [
    "#7DD06F",
    "#844081",
    "#688EC1",
    "#C17E73",
    "#484125",
    "#6CD3A7",
    "#597873",
    "#7B6FD0",
    "#CF4A31",
    "#D0CD47",
    "#722A2D",
    "#CBC594",
    "#D19EC4",
    "#5A7E36",
    "#D4477D",
    "#403552",
    "#76D73C",
    "#96CED5",
    "#CE54D1",
    "#C48736",
]

kelly: List[str] = [
    "#FFB300",
    "#803E75",
    "#FF6800",
    "#A6BDD7",
    "#C10020",
    "#CEA262",
    "#817066",
    "#007D34",
    "#F6768E",
    "#00538A",
    "#FF7A5C",
    "#53377A",
    "#FF8E00",
    "#B32851",
    "#F4C800",
    "#7F180D",
    "#93AA00",
    "#593315",
    "#F13A13",
    "#232C16",
]
bear: List[str] = [
    "#faa818",
    "#41a30d",
    "#fbdf72",
    "#367d7d",
    "#d33502",
    "#6ebcbc",
    "#37526d",
    "#916848",
    "#f5b390",
    "#342739",
    "#bed678",
    "#a6d9ee",
    "#0d74b6",
    "#60824f",
    "#725ca5",
    "#e0598b",
]

ironMan: List[str] = [
    "#371377",
    "#7700FF",
    "#9E0142",
    "#FF0080",
    "#DC494C",
    "#F88D51",
    "#FAD510",
    "#FFFF5F",
    "#88CFA4",
    "#238B45",
    "#02401B",
    "#0AD7D3",
    "#046C9A",
    "#A2A475",
    "grey35",
]

circus: List[str] = [
    "#D52126",
    "#88CCEE",
    "#FEE52C",
    "#117733",
    "#CC61B0",
    "#99C945",
    "#2F8AC4",
    "#332288",
    "#E68316",
    "#661101",
    "#F97B72",
    "#DDCC77",
    "#11A579",
    "#89288F",
    "#E73F74",
]

paired: List[str] = [
    "#A6CDE2",
    "#1E78B4",
    "#74C476",
    "#34A047",
    "#F59899",
    "#E11E26",
    "#FCBF6E",
    "#F47E1F",
    "#CAB2D6",
    "#6A3E98",
    "#FAF39B",
    "#B15928",
]

grove: List[str] = [
    "#1a1334",
    "#01545a",
    "#017351",
    "#03c383",
    "#aad962",
    "#fbbf45",
    "#ef6a32",
    "#ed0345",
    "#a12a5e",
    "#710162",
    "#3B9AB2",
]

summerNight: List[str] = [
    "#2a7185",
    "#a64027",
    "#fbdf72",
    "#60824f",
    "#9cdff0",
    "#022336",
    "#725ca5",
]

zissou: List[str] = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]
darjeeling: List[str] = ["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"]
rushmore: List[str] = ["#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F"]
captain: List[str] = ["grey", "#A1CDE1", "#12477C", "#EC9274", "#67001E"]
horizon: List[str] = [
    "#000075",
    "#2E00FF",
    "#9408F7",
    "#C729D6",
    "#FA4AB5",
    "#FF6A95",
    "#FF8B74",
    "#FFAC53",
    "#FFCD32",
    "#FFFF60",
]
horizonExtra: List[str] = [
    "#000436",
    "#021EA9",
    "#1632FB",
    "#6E34FC",
    "#C732D5",
    "#FD619D",
    "#FF9965",
    "#FFD32B",
    "#FFFC5A",
]
blueYellow: List[str] = [
    "#352A86",
    "#343DAE",
    "#0262E0",
    "#1389D2",
    "#2DB7A3",
    "#A5BE6A",
    "#F8BA43",
    "#F6DA23",
    "#F8FA0D",
]
sambaNight: List[str] = [
    "#1873CC",
    "#1798E5",
    "#00BFFF",
    "#4AC596",
    "#00CC00",
    "#A2E700",
    "#FFFF00",
    "#FFD200",
    "#FFA500",
]
solarExtra: List[str] = [
    "#3361A5",
    "#248AF3",
    "#14B3FF",
    "#88CEEF",
    "#C1D5DC",
    "#EAD397",
    "#FDB31A",
    "#E42A2A",
    "#A31D1D",
]
whitePurple: List[str] = [
    "#f7fcfd",
    "#e0ecf4",
    "#bfd3e6",
    "#9ebcda",
    "#8c96c6",
    "#8c6bb1",
    "#88419d",
    "#810f7c",
    "#4d004b",
]
whiteBlue: List[str] = [
    "#fff7fb",
    "#ece7f2",
    "#d0d1e6",
    "#a6bddb",
    "#74a9cf",
    "#3690c0",
    "#0570b0",
    "#045a8d",
    "#023858",
]
whiteRed: List[str] = ["white", "red"]
comet: List[str] = ["#E6E7E8", "#3A97FF", "#8816A7", "black"]

greenBlue: List[str] = [
    "#e0f3db",
    "#ccebc5",
    "#a8ddb5",
    "#4eb3d3",
    "#2b8cbe",
    "#0868ac",
    "#084081",
]

beach: List[str] = ["#87D2DB", "#5BB1CB", "#4F66AF", "#F15F30", "#F7962E", "#FCEE2B"]
coolwarm: List[str] = ["#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"]
fireworks: List[str] = ["white", "#2488F0", "#7F3F98", "#E22929", "#FCB31A"]
greyMagma: List[str] = ["grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"]
fireworks2: List[str] = ["black", "#2488F0", "#7F3F98", "#E22929", "#FCB31A"]
purpleOrange: List[str] = ["#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"]
