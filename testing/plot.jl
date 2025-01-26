using Pkg

Pkg.add(["GLMakie", "CSV", "DataFrames"])

using GLMakie
using CSV
using DataFrames

function plot_quaternion_data(filename)
    # Read the CSV file
    df = CSV.read(filename, DataFrame)

    # Create the figure and axis
    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1],
        xlabel="X",
        ylabel="Y",
        zlabel="Z",
        title="3D Visualization of Quaternion Data"
    )

    # Scatter plot with color mapping
    scatter!(ax,
        df.x, df.y, df.z,
        color=df.scalar,
        colormap=:viridis,
        markersize=0.1
    )

    # Set axis limits
    #  xlims!(ax, -1000, 1000)
    #  ylims!(ax, -1000, 1000)
    #  zlims!(ax, 0, 1000)

    # Add colorbar
    Colorbar(fig[1, 2], colormap=:viridis)

    # Display the figure
    display(fig)
end

# Usage
plot_quaternion_data("output.csv")
