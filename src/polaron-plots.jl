using CairoMakie
using LaTeXStrings

include("data-processing.jl")

function add_labels(df)
    os = df.is_optimized
    Ns = df.N
    gs = df.g
    labels = LaTeXString[]
    if hasproperty(df, :u)
        us = df.u
    else
        us = zeros(size(df, 1))
    end

    n_os = length(unique(os))
    n_Ns = length(unique(Ns))
    n_us = length(unique(us))

    labels = map(zip(os, Ns, gs, us)) do (o, N, g, u)
        m = round(Int, √N)
        label = "\$$m×$m"
        if n_us > 1
            label = string(label, ",\\,U_{IB}=$u")
        end
        if n_os > 1
            if o
                label = ",\\,(g\$ opt.\$)"
            else
                label = string(label, ",\\,(g=$g)")
            end
        end
        if N ≠ m^2
            label = string(label, ",\\,N=$N")
        end
        LaTeXString(string(label, "\$"))
    end
    df[!,:label] = labels
    return df
end

function plot_ep(df; pe=false, showsus=false, debug=nothing, debug_name=string(debug), kwargs...)
    if !hasproperty(df, :Ep_S)
        df = process_polaron(df)
    end
    df = add_labels(df)
    fig = Figure()
    ax1 = Axis(fig[1, 1]; xlabel=L"4t/U", ylabel=L"(E_{P} + 4t)/U", kwargs...)
    if !isnothing(debug)
        ax2 = Axis(fig[2, 1]; xlabel=L"4t/U", ylabel=debug_name)
    end

    groups = groupby(df, [:N, :is_optimized])
    sort!(df, [:N])

    for (i, d) in enumerate(groups)
        label = d.label[1]
        marker = [s ? :star : :circle for s in d.is_suspicious]
        if !pe
            errorbars!(ax1, 4 * d.t, d.Ep_S, d.Ep_S_err; label, marker)
            if showsus
                ok = @rsubset d !:is_suspicious
                sus = @rsubset d :is_suspicious
                scatter!(ax1, 4 * sus.t, sus.Ep_S; label, marker=:utriangle, color=Cycled(i))
            else
                ok = d
            end
            scatter!(ax1, 4 * ok.t, ok.Ep_S; label, color=Cycled(i))
        else
            errorbars!(ax1, 4 * d.t, d.Ep_pe, d.Ep_pe_err_l, d.Ep_pe_err_u; label)
            scatter!(ax1, 4 * d.t, d.Ep_pe; label)
        end
        if !isnothing(debug)
            scatter!(ax2, 4 * d.t, d[:,debug]; label, color=Cycled(i))
        end
    end

    if length(groups) > 1
        Legend(fig[1, 2], ax1; merge=true)
    end

    return fig
end

function plot_debug(df)
    if hasproperty(df, :Ep)
        error("wrong df")
    end
    df = add_labels(df)
    fig = Figure(resolution=(800,800))
    ax1 = Axis(fig[1:2, 1:2]; xlabel=L"4t/U", ylabel=L"⟨S⟩ / U")
    ax2 = Axis(fig[3:4, 1:2]; xlabel=L"4t/U", ylabel=L"E_p / U")
    ax3 = Axis(fig[5, 1]; xlabel=L"4t/U", ylabel=L"pcb$$")
    ax4 = Axis(fig[5, 2]; xlabel=L"4t/U", ylabel=L"Ω")
    ax5 = Axis(fig[6, 1]; xlabel=L"4t/U", ylabel=L"ref. occ.$$")
    ax6 = Axis(fig[6, 2]; xlabel=L"4t/U", ylabel=L"len / norm$$")

    groups = groupby(df, [:N, :is_optimized, :u])
    markers = [
        :circle, :rect, :diamond, :utriangle, :dtriangle, :ltriangle, :rtriangle
    ]

    vlines!(ax1, [4 * t_crit]; color=Cycled(2), linestyle=:dot)

    for (i, d) in enumerate(groups)
        label = d.label[1]
        marker = markers[mod1(i, length(markers))]
        errorbars!(ax1, 4 * d.t, d.S_mean, d.S_err; label)
        scatter!(ax1, 4 * d.t, d.S_mean; label, marker)

        errorbars!(ax2, 4 * d.t, d.pe_val, d.pe_err_l, d.pe_err_u; label)
        scatter!(ax2, 4 * d.t, d.pe_val; label, marker)

        errorbars!(ax3, 4 * d.t, d.bias, d.bias_err; label)
        scatter!(ax3, 4 * d.t, d.bias; label, marker)

        scatter!(ax4, 4 * d.t, d.Ω; label, marker)

        scatter!(ax5, 4 * d.t, d.reference_occupation; label, marker)
        scatter!(ax6, 4 * d.t, d.len_vs_norm; label, marker)
    end

    Legend(fig[:, 2], ax1; merge=true)

    return fig
end

function plot_series(dfs...; key2=nothing, window=1000)
    fig = Figure()
    ax1 = Axis(fig[1, 1]; xlabel=L"step$$", ylabel=L"S")
    if !isnothing(key2)
        ax2 = Axis(fig[2,1]; xlabel=L"step$$", ylabel=string(key2))
    end
    for d in dfs
        if key2 == :cum_mean
            d[!,:cum_mean] = reverse([mean(d.shift[i:end]) for i in eachindex(d.shift)])
        end
        label = LaTeXString(string(
            "\$N=", get(metadata(d), "N", ""),
            ", t=", get(metadata(d), "t", ""),
            ", U_{IB}=", get(metadata(d), "u", ""),
            "\$",
        ))
        if key2 == :sliding_mean
            d[!, :sliding_mean] = [
                mean(d.shift[max(i-50000,1):i]) for i in eachindex(d.shift)
                    ]
        end

        lines!(ax1, d.steps, d.shift; label)
        if !isnothing(key2)
            lines!(ax2, d.steps, d[:,key2]; label)
        end
    end
    Legend(fig[:,2], ax1)
    return fig
end

function sliding_mean(dfs...; reverse=false, step=100, skip=10)
    if reverse
        title=L"Skipping from the front$$"
    else
        title=L"Skipping from the back$$"
    end

    fig = Figure()
    ax1 = Axis(fig[1,1]; xlabel=L"step$$", ylabel=L"S - μS_{end}", title)
    ax2 = Axis(fig[2,1]; xlabel=L"Ω", ylabel=L"⟨S⟩ - μS_{end}")
    ax3 = Axis(fig[3,1]; xlabel=L"Ω", ylabel=L"PCB$$")

    for d in dfs
        label = LaTeXString(string(
            "\$N=", get(metadata(d), "N", ""),
            ", t=", get(metadata(d), "t", ""),
            ", U_{IB}=", get(metadata(d), "u", ""),
            "\$",
        ))

        pd = DataFrame()
        S_end = mean(d.shift)
        for i in 0:step:size(d, 1)
            if reverse
                curr_shift = d.shift[(i + 1):end]
            else
                curr_shift = d.shift[1:(end-i)]
            end
            b = blocking_analysis(curr_shift)
            S_mean = b.mean .- S_end
            S_err = b.err
            S_err_err = b.err_err
            bias = S_err^2/2 * sum(df.dτ)
            bias_err = S_err_err * S_err * sum(df.dτ)

            push!(pd, (
                ; Ω=length(curr_shift), S_mean, S_err, bias, bias_err
            ))
        end

        pd = sort!(pd, [:Ω])[(skip+1):end,:]

        lines!(ax1, d.steps, d.shift .- S_end; label)
        band!(ax2, pd.Ω, pd.S_mean - pd.S_err, pd.S_mean + pd.S_err; label)
        lines!(ax2, pd.Ω, pd.S_mean; label)
        band!(ax3, pd.Ω, pd.bias - pd.bias_err, pd.bias + pd.bias_err; label)
        lines!(ax3, pd.Ω, pd.bias; label)
    end

    hlines!(ax1, [0]; label=L"y=0", color=:gray)
    Legend(fig[:,2], ax1; merge=true)

    return fig
end

function segmented_mean(dfs...; step=1000, data=true)
    fig = Figure()
    ax1 = Axis(fig[1, 1])

    i = 1
    for d in dfs
        label = LaTeXString(string(
            "\$N=", get(metadata(d), "N", ""),
            ", t=", get(metadata(d), "t", ""),
            ", U_{IB}=", get(metadata(d), "u", ""),
            "\$",
        ))
        if data
            lines!(ax1, d.steps, d.shift; label, color=Cycled(i))
        end

        xs = Int[]
        ys = Float64[]
        es = Float64[]

        for i in 0:step:(size(d, 1) - step)
            shift = d.shift[i .+ (1:step)]
            b = blocking_analysis(shift)
            append!(xs, [i+1, i+step, i+step])
            append!(ys, [b.mean, b.mean, NaN])
            append!(es, [b.err, b.err, NaN])
        end

        lines!(ax1, xs, ys .- es; label, color=Cycled(i+1), linestyle=:dot)
        lines!(ax1, xs, ys .+ es; label, color=Cycled(i+1), linestyle=:dot)
        lines!(ax1, xs, ys; label, color=Cycled(i+1))
        i += 2
    end

    return fig
end
