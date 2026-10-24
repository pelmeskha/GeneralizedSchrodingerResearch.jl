function live_plot_record(run_iterations, x, U0; filename="output.mp4", framerate=30)
    observable_U = Observable(U0)
    observable_x = Observable(x)
    figure = Figure()
    axis = Axis(figure[1, 1], limits = (minimum(x), maximum(x), 0.0, 1.0))
    lines!(axis, observable_x, observable_U)
    display(figure)
    record(figure, filename, framerate=framerate) do io
        show_frame(t, U) = (observable_U[] = U; recordframe!(io))
        run_iterations(show_frame)
    end
end
