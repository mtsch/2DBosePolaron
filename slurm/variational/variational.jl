using Gutzwiller
using Rimu
using EasyArgs

function setup_H(m; N=m*m, t=0.1, u_ib=0.0, g=0.0)
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

function main(
    ;
    m=get_arg("m", 6),
    N=get_arg("N", m * m),
    t=get_arg("t", 0.1),
    u_ib=get_arg("u_ib", 0.2),
    steps=get_arg("steps", 1e6),
    warmup=get_arg("steps", 1e6),
    tasks=get_arg("tasks", 100),
    )

    @info "Starting" m N t u_ib steps warmup tasks

    H = setup_H(m; t, u_ib)
    res = gutz_optimize(
        H, 1.0;
        step=0.8, g_tol=0.01, qmc=true, verbose=true,
        warmup, steps, tasks,
    )
    println(stdout, join((
        m, N, t, u_ib, warmup, steps, tasks,
        res.minimizer, res.step, res.minimum.mean, res.minimum.err,
        res.minimum.acceptance, res.converged, res.evals, res.iters,
    ), ","))
end

if !isinteractive()
    main()
end
