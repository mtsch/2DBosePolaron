using Rimu
using EasyArgs
using ProgressMeter

function main(sz=get_arg("s", 200_000))
    @showprogress for f in readdir(".")
        if endswith(f, ".arrow")
            df = load_df(f)
            if size(df, 1) == sz
                @info "rewriting `$f`"
                save_df(f, df)
            end
        end
    end
    @info "done."
end

if !isinteractive()
    main()
end
