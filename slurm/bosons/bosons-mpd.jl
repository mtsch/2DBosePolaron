using Rimu
using EasyArgs
using DataFrames
using Rimu.RMPI

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
    return DVec(zip(addrs, ones(length(addrs))); style=IsDynamicSemistochastic())
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
        H = HubbardRealSpace(addr; geometry, u=1.0, t=1.0)
    end

    if g ≠ 0
        H = GutzwillerSampling(H, g)
    end
    return H
end

function save_vec(f, v)
    ks = collect(keys(v))
    vs = collect(values(v))
    RimuIO.save_df(f, DataFrame(keys=ks, values=vs))
end

function load_vec(f)
    df = RimuIO.load_df(f)
    return DVec(zip(df.keys, df.values); style=IsDynamicSemistochastic())
end

function warmup(H, v_init, steps_warmup, Nt, dτ, v_filename)
    @mpi_root @info "Warming up..."
    el = @elapsed begin
        norm1 = walkernumber(v_init)
        if norm1 > Nt
            v_init = map!(v_init, values(v_init)) do x
                x * Nt / norm1
            end
            shift_init = rayleigh_quotient(H, v_init)
        else
            shift_init = 0.0
        end

        params = RunTillLastStep(; shift=shift_init, dτ, laststep=steps_warmup)
        s_strat = DoubleLogUpdateAfterTargetWalkers(targetwalkers=Nt)
        r_strat = ReportDFAndInfo(; writeinfo=true, io=stderr, info_interval=1000)

        if v_filename == ""
            v = MPIData(copy(v_init))
        else
            el_load = @elapsed v = load_dvec(v_filename)
            @info "Vector loaded from $v_filename in $el_load seconds..." walkernumber(v)
        end

        df, st = lomc!(H, v; s_strat, params, r_strat)
        while params.step ≠ params.laststep
            dτ /= 10
            @mpi_root @info "setting dτ to $dτ"
            v = MPIData(copy(v_init))
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
    #occ_counter = map!(similar(v_init), values(v_init)) do x
    #    1 / len
    #end
    post_step = (
        ProjectedEnergy(H, normalize(v_init)),
        Projector(norm2=Norm2Projector()),
     #   Projector(reference_occupation=occ_counter),
        Rimu.Timer(),
    )

    params = st.replicas[1].params
    params.step = 0
    params.laststep = steps_measure
    v = st.replicas[1].v
    filename *= ".arrow"
    r_strat = ReportToFile(; filename, chunk_size, io=stderr)
    s_strat = DoubleLogUpdateAfterTargetWalkers(targetwalkers=Nt)

    _, st = lomc!(H, v; params, post_step, r_strat, s_strat)
    if v_filename ≠ ""
        @info "Saving vector..."
        v = st.replicas[1].v
        el = @elapsed save_dvec(v_filename, v)
        @info "Done in $el seconds"
    end
    @mpi_root @info "Done"
end

function main(
    ;
    m=get_arg("m", 6),
    N=m*m,
    t=get_arg("t", 0.1),
    u_ib=get_arg("u", 0.2),
    g=get_arg("g", 0.2),
    Nt=Int(get_arg("N", 1e7)),
    dτ=get_arg("dt", 0.0001),
    steps_warmup=get_arg("warmup", 25_000),
    steps_measure=get_arg("measure", 25_000),
    filename=get_arg("o", "test"),
    chunk_size=get_arg("chunk_size", 1000),
    save_dvec_to=get_arg("save", ""),
    initial_dvec=get_arg("i", ""),
    )
    EasyArgs.check_unused_args(; error=true)
    @mpi_root @info "Starting. M=$m × $m, N=$N" t u_ib g Nt dτ steps_warmup steps_measure filename mpi_size() save_dvec_to initial_dvec

    id = get(ENV, "SLURM_JOB_ID", "")
    filename = string(filename, "_N", N, "_m", m, "_t", t, "_u", u_ib, "_Nt", Nt, "_g", g, "_", id)

    @mpi_root begin
        meta = filename * ".txt"
        open(meta, "w") do f
            println(f, "filename,m,N,t,u_ib,g,Nt,dτ,steps_warmup,steps_measure")
            println(f, join((filename,m,N,t,u_ib,g,Nt,dτ,steps_warmup,steps_measure),","))
        end
    end

    H = setup_H(m, N, t, u_ib, g)
    v_init = setup_vec(m, N, u_ib≠0)
    st = warmup(H, v_init, steps_warmup, Nt, dτ, initial_dvec)
    measure(H, v_init, st, Nt, steps_measure, filename, chunk_size, save_dvec_to)
end

if !isinteractive()
    mpi_allprintln("hello")
    main()
end
