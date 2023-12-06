include("../slurm/bosons/bosons.jl")
using Rimu
using DataFramesMeta
using KrylovKit
using ProgressMeter

function mini_energy(m, t, u_ib, ne)
    ham_imp = setup_H(m, m*m, t, u_ib, 0)
    sp = sparse(ham_imp)
    eigsolve(sp, ne, :SR)[1][1:ne]
end

function mini_series(m=3, ne=1; ts=0.01:0.001:0.2)
    df = DataFrame()
    @showprogress for t in ts
        for u_ib in (0.2)
            E_noi = mini_energy(m, t, 0.0, ne)
            E_imp = mini_energy(m, t, u_ib, ne)

            for i in 1:ne
                push!(df, (; m, t, u_ib, i=i-1, E_noi=E_noi[i], E_imp=E_imp[i], E_p=E_imp[i]-E_noi[i]))
            end
        end
    end
    return df
end

function mini2ep(df)
    @chain df begin
        @rsubset :u_ib == 0.2
        @rtransform :N = :m * :m :g = 0.0 :is_optimized=false :Ep_err = 0.0 :Ep = :E_p + 4*:t
        @select :N :t :Ep :Ep_err
    end
end
