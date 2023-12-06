using Rimu
using Gutzwiller
using EasyArgs
using DataFrames
using Rimu.RMPI
using Arrow

function setup_vec(m, N, impurity=true)
    M = m*m
    if impurity
        addrs = [
            CompositeFS(
                near_uniform(BoseFS{N,M}),
                BoseFS(M, i => 1)
            ) for i in 1:M
        ]
    else
        addrs = [near_uniform(BoseFS{N,M})]
    end
    return PDVec(zip(addrs, ones(length(addrs))); style=IsDynamicSemistochastic())
end

function setup_H(m, N, t, u_ib, g)
    M = m * m
    geometry = PeriodicBoundaries(m, m)
    if u_ib ≠ 0
        u = [1 u_ib; u_ib 0]
        t = [t, t]
        addr = CompositeFS(
            near_uniform(BoseFS{N,M}),
            BoseFS(M, 1 => 1),
        )
        H = HubbardRealSpace(addr; geometry, u, t)
    else
        addr = near_uniform(BoseFS{N,M})
        H = HubbardRealSpace(addr; geometry, u=[1.0], t=[t])
    end

    if g ≠ 0
        H = GutzwillerSampling(H, g)
    end
    return H
end

function save_state(f, v, shift, dτ)
    metadata = ("S" => string(shift), "dt" => string(dτ))
    Arrow.write(
        f, (; keys=collect(keys(v)), values=collect(values(v)));
        compress=:zstd, metadata,
    )
end

function load_state(f)
    tbl = Arrow.Table(f)
    shift = parse(Float64, Arrow.metadata(tbl)[]["S"])
    dτ = parse(Float64, Arrow.metadata(tbl)[]["dt"])
    return PDVec(zip(tbl.keys, tbl.values); style=IsDynamicSemistochastic()), shift, dτ
end

function warmup(H, v_init, steps_warmup, Nt, dτ, v_filename)
    @mpi_root @info "Warming up..."
    el = @elapsed begin
        norm1 = walkernumber(v_init)

        if v_filename == ""
            v = copy(v_init)
            if norm1 > Nt
                v_init = map!(v_init, values(v_init)) do x
                    x * Nt / norm1
                end
                shift_init = rayleigh_quotient(H, v_init)
            else
                shift_init = 0.0
            end
            s_strat = DoubleLogUpdate(targetwalkers=Nt)
        else
            el_load = @elapsed v, shift_init, dτ = load_state(v_filename)
            dτ = dτ * 10
            @info "Vector loaded from $v_filename in $el_load seconds." walkernumber(v) shift_init dτ
            s_strat = DoubleLogUpdate(targetwalkers=Nt)
        end

        params = RunTillLastStep(; shift=shift_init, dτ, laststep=steps_warmup)
        r_strat = ReportDFAndInfo(; writeinfo=is_mpi_root(), io=stderr, info_interval=1000)

        df, st = lomc!(H, v; s_strat, params, r_strat)
        while params.step ≠ params.laststep
            dτ /= 2
            @mpi_root @info "setting dτ to $dτ"
            v = copy(v_init)
            params = RunTillLastStep(; shift=shift_init, dτ, laststep=steps_warmup)
            df, st = lomc!(H, v; params, s_strat, r_strat)
        end
    end
    v = st.replicas[1].v
    wn_final = walkernumber(v)
    @mpi_root @info "Took $el seconds with a final norm of $(wn_final)"
    return st
end

function measure(H, v_init, st, Nt, steps_measure, filename, chunk_size, v_filename)
    @mpi_root @info "Running..."
    len = length(v_init)
    occ_counter = map!(similar(v_init), values(v_init)) do x
        1 / len
    end
    post_step = (
        ProjectedEnergy(H, normalize(v_init)),
        Projector(norm2=Norm2Projector()),
        Projector(reference_occupation=occ_counter),
        Rimu.Timer(),
    )

    params = st.replicas[1].params
    params.step = 0
    params.laststep = steps_measure
    v = st.replicas[1].v
    filename *= ".arrow"
    r_strat = ReportToFile(; filename, chunk_size, io=stderr)
    s_strat = DoubleLogUpdate(targetwalkers=Nt)

    _, st = lomc!(H, v; params, post_step, r_strat, s_strat)
    if v_filename ≠ ""
        @info "Saving vector..."
        v = st.replicas[1].v
        shift = st.replicas[1].params.shift
        dτ = st.replicas[1].params.dτ
        el = @elapsed save_state(v_filename, v, shift, dτ)
        @info "Done in $el seconds"
    end

    @mpi_root begin
        @info "Resaving file"
        df = RimuIO.load_df(r_strat.filename)
        RimuIO.save_df(r_strat.filename, df)
        @info "Done"
    end
end

function main(
    ;
    m=get_arg("m", 6),
    N=get_arg("nparts", m*m),
    t=get_arg("t", 0.1),
    u_ib=get_arg("u", 0.2),
    Nt=Int(get_arg("N", 1e7)),
    dτ=get_arg("dt", 0.512),
    steps_warmup=get_arg("warmup", 25_000),
    steps_measure=get_arg("measure", 25_000),
    filename=get_arg("o", "test"),
    chunk_size=get_arg("chunk_size", 1000),
    save_dvec_to=get_arg("save", ""),
    initial_dvec=get_arg("load", ""),
    )
    EasyArgs.check_unused_args(; error=true)
    @mpi_root @info "Starting. M=$m × $m, N=$N" t u_ib Nt dτ steps_warmup steps_measure filename mpi_size() save_dvec_to initial_dvec

    id = get(ENV, "SLURM_JOB_ID", "")

    @mpi_root @info "Optimizing g"
    N_small = round(Int, 9N/(m*m))
    H_base = setup_H(3, round(Int, N_small), t, u_ib, 0.0)
    res = gutz_optimize(H_base, 1.0; verbose=false)
    g = res.minimizer
    @mpi_root @info "Done." g

    filename = string(filename, "_optimized_N", N, "_m", m, "_t", t, "_u", u_ib, "_Nt", Nt, "_g", g, "_", id)

    H = setup_H(m, N, t, u_ib, g)
    v_init = setup_vec(m, N, u_ib≠0)
    st = warmup(H, v_init, steps_warmup, Nt, dτ, initial_dvec)
    measure(H, v_init, st, Nt, steps_measure, filename, chunk_size, save_dvec_to)
end

if !isinteractive()
    mpi_allprintln("hello")
    main()
end
