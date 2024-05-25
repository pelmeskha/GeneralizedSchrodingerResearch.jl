function set_grid_alpha(figure, value::Float64)
    0.0 ≤ value ≤ 1.0 || throw(ArgumentError(":gridalpha parameter must be in [0,1]"))
    xaxis = Plots.get_axis(Plots.get_subplot(figure,1),:x)
    yaxis = Plots.get_axis(Plots.get_subplot(figure,1),:y)
    xaxis[:gridalpha] = value
    yaxis[:gridalpha] = value
end
