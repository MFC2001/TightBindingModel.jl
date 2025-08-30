
function _plot_init()
    
    # 检查Plots是否可用
    if isdefined(Main, :Plots)
        include("./Plots/Plots.jl")
        @info "Loaded Plots support."
    end
    if isdefined(Main, :CairoMakie)
        include("./CairoMakie/CairoMakie.jl")
        @info "Loaded CairoMakie support."
    end
    if isdefined(Main, :GLMakie)
        include("./GLMakie/GLMakie.jl")
        @info "Loaded GLMakie support."
    end
    # 检查轻量级替代方案
    if isdefined(Main, :UnicodePlots)
        include("./UnicodePlots/UnicodePlots.jl")
        @info "Using UnicodePlots support."
    end

end

_plot_init()