using Rimu
using DataFramesMeta

const t_crit = 5.909e-2

function load_data(f::String)
    filename = joinpath(@__DIR__, "../data", f)
    df = RimuIO.load_df(filename)
    if !hasproperty(df, :time_h)
        @info "re-saving `$f`"
        df = @rtransform df :time_h = (:time - df.time[1]) / 3600
        RimuIO.save_df(filename, df)
    end
    for (k, v) in pairs(filename_metadata(f))
        metadata!(df, string(k), v)
    end
    return df
end

function load_data(; datadir=joinpath(@__DIR__, "../data"), kwargs...)
    all_files = readdir(datadir)

    matches = filter(all_files) do f
        meta = filename_metadata(f)
        !isnothing(meta) && all(kwargs) do (k, v)
            meta[k] == v
        end
    end
    if length(matches) == 1
        load_data(only(matches))
    elseif length(matches) == 0
        @warn "No matches found"
        return nothing
    else
        @warn "Multiple matches found" matches
        return nothing
    end
end

function convergence_test(shift, id)
    μ1, σ1 = mean_and_se(shift[1:(end ÷ 2)])
    μ2, σ2 = mean_and_se(shift[(end ÷ 2 + 1):end])

    if abs(μ1 - μ2) > σ1 + σ2
        @warn "data with id $id is suspicious" abs(μ1 - μ2) 3σ1 3σ2
        return true
    else
        return false
    end
end

function filename_metadata(f)
    id_match = match(r"([0-9]+)\.arrow", f)
    if !isnothing(id_match)
        id = id_match[1]
    else
        id = missing
    end
    try
        g = parse(Float64, match(r"_g([0-9.]+)", f)[1])

        return (
            ;
            g = parse(Float64, match(r"_g([0-9.]+)", f)[1]),
            N = parse(Int, match(r"_N([0-9]+)", f)[1]),
            t = parse(Float64, match(r"_t([0-9.]+)", f)[1]),
            u = parse(Float64, match(r"_u([0-9.]+)", f)[1]),
            Nt = parse(Float64, match(r"_Nt([0-9]+)", f)[1]),
            is_optimized = match(r"optimized", f) ≢ nothing,
            id,
        )
    catch e
        return nothing
    end
end

function process_group(name; datadir=joinpath(@__DIR__, "../data"), skip=0, force=false)
    res = DataFrame()
    filename = joinpath(datadir, string(name, ".arrow"))
    if !force && isfile(filename)
        @info "Loading from $filename."
        return RimuIO.load_df(filename)
    elseif force
        @info "Forcing update..."
    end

    el = @elapsed for f in readdir(datadir)
        if startswith(f, name) && endswith(f, ".arrow") && joinpath(datadir, f) ≠ filename
            println(f)
            meta = filename_metadata(f)
            df = load_df(joinpath(datadir, f))[(skip + 1):end,:]
            if size(df, 1) == 0
                continue
            end

            is_suspicious = convergence_test(df.shift, meta.id)
            se = shift_estimator(df)
            S_mean = se.mean
            S_err = se.err
            S_err_err = se.err_err

            pe_val, pe_err_l, pe_err_u = try
                pe = projected_energy(df)
                pes = val_and_errs(pe)
                if !isfinite(pes[1])
                    NaN, NaN, NaN
                else
                    pes
                end
            catch e
                NaN, NaN, NaN
            end

            bias = S_err^2/2 * sum(df.dτ)
            bias_err = S_err_err * S_err * sum(df.dτ)

            reference_occupation = mean(df.reference_occupation)
            len_vs_norm = mean(df.len ./ df.norm)

            push!(res,
                  (;
                   N=meta.N, g=meta.g, is_optimized=meta.is_optimized,
                   u=meta.u, t=meta.t, Nt=meta.Nt, Ω=length(df.shift), dτ=df.dτ[1],
                   S_mean, S_err, pe_val, pe_err_l, pe_err_u,
                   bias, bias_err,
                   reference_occupation, len_vs_norm, filename=f,
                   is_suspicious,
                   ); promote=true)
        end
    end
    sort!(res, [:N, :is_optimized, :u, :t])
    @info "Saving to $filename."
    RimuIO.save_df(filename, res)
    @info "Took $el seconds."
    return res
end

function process_polaron(df)
    df_noi = @chain df begin
        @rsubset :u == 0
        @select begin
            :N
            :is_optimized
            :g
            :t
            :S_mean_noi = :S_mean
            :S_err_noi = :S_err
            :pe_val_noi = :pe_val
            :pe_err_u_noi = :pe_err_u
            :pe_err_l_noi = :pe_err_l
            :Ω_noi = :Ω
            :len_vs_norm_noi = :len_vs_norm
            :bias_noi = :bias
            :bias_err_noi = :bias_err
            :is_suspicious_noi = :is_suspicious
        end
    end
    df_imp = @chain df begin
        @rsubset :u == 0.2
        @select begin
            :N
            :is_optimized
            :t
            :S_mean_imp = :S_mean
            :S_err_imp = :S_err
            :pe_val_imp = :pe_val
            :pe_err_u_imp = :pe_err_u
            :pe_err_l_imp = :pe_err_l
            :Ω_imp = :Ω
            :len_vs_norm_imp = :len_vs_norm
            :bias_imp = :bias
            :bias_err_imp = :bias
            :is_suspicious_imp = :is_suspicious
        end
    end

    df_proc = @chain innerjoin(df_noi, df_imp, on=[:N, :is_optimized, :t]) begin
        @rtransform begin
            :Ep_S = (:S_mean_imp - :S_mean_noi) + 4*:t
            :Ep_S_err = √(:S_err_noi^2 + :S_err_imp^2)
            :Ω = min(:Ω_imp, :Ω_noi)
            :len_vs_norm = max(:len_vs_norm_imp, :len_vs_norm_noi)
            :Ep_pe = (:pe_val_imp - :pe_val_noi) + 4*:t
            :Ep_pe_err_u = √(:pe_err_u_noi^2 + :pe_err_u_imp^2)
            :Ep_pe_err_l = √(:pe_err_l_noi^2 + :pe_err_l_imp^2)
            :bias = :bias_noi + :bias_imp
            :bias_err = √(:bias_err_noi^2 + :bias_err_imp^2)
            :is_suspicious = :is_suspicious_imp || :is_suspicious_noi
        end
        @select begin
            :N
            :g
            :is_optimized
            :t
            :Ep_S
            :Ep_S_err
            :Ep_pe
            :Ep_pe_err_u
            :Ep_pe_err_l
            :Ω
            :len_vs_norm
            :bias
            :bias_err
            :is_suspicious
        end
    end
    return df_proc
end
