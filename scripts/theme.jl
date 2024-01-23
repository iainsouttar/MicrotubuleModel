Attributes(
    resolution=(1500,600),
    fontsize=30,
    font="monserrat",
    palette = (
        color = colorschemes[:seaborn_colorblind].colors,
        marker = [:rect, :circle, :utriangle, :diamond], 
        linestyle = [:dash, :dot, :dashdot, :dashdotdot],
    ),
    linewidth=3,
    marker='◯',
    markersize=25,
    strokewidth=1.5,
    strokecolor="#DD8452",
    colormap=:viridis,
    Lines = (color = "#4C72B0", linewidth=3),
    Scatter = (color="#55A868", marker='◯', markersize=15, strokewidth=1.5, strokecolor="#55A868"),
    grid=false,
    Axis = (
        xgridvisible = false,
        ygridvisible = false,
    )
)