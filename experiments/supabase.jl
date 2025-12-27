using HTTP
using JSON
using CSV
using DataFrames
using DeRham

include("find_newton_polygons/cpu_experiments.jl")

############################################################
# SUPABASE CODE
############################################################

struct SupabaseClient
    url::String
    key::String
end

function headers(c::SupabaseClient)
    return [
        "apikey"        => c.key,
        "Authorization" => "Bearer $(c.key)",
        "Content-Type"  => "application/json"
    ]
end

function push_row(c::SupabaseClient, table::String, row::Dict)
    baseurl = "$(c.url)/rest/v1/$table"
    payload = JSON.json(row)

    HTTP.post(
        baseurl,
        headers = [
            "apikey"        => c.key,
            "Authorization" => "Bearer $(c.key)",
            "Content-Type"  => "application/json"
        ],
        body = payload
    )
end

function fetch_table(c::SupabaseClient, table::String)

    resp = nothing
    try
        url = "$(c.url)/rest/v1/$table?select=*"
        resp = HTTP.get(url, headers=headers(c))
    catch e
        println("Error fetching table: $(e)")
        error("Please make the table $(table) first on Supabase.")
    end

    rows = JSON.parse(String(resp.body))

    if length(rows) == 0
        return DataFrame(
            n = Int[],
            d = Int[],
            p = Int[],
            slopes = String[],
            slopelengths = String[],
            values = String[],
            slopesbefore = String[],
            polystr = String[]
        )
    end
    return DataFrame(rows)
end

############################################################
# ENCODERS
############################################################

encode_key(np) = JSON.json(np)
encode_value(f) = string(f)

function update_df(df::DataFrame, n, d, p, np, f, push_to_supabase=false)
    
    if nrow(df) > 0
        existing = df[
            (df.n .== n) .& 
            (df.d .== d) .& 
            (df.p .== p) .& 
            (df.values .== string(np.values)), :]
        if nrow(existing) > 0
            return df
        end
    end
    
    row = Dict(
        "n"            => n,
        "d"            => d,
        "p"            => p,
        "slopes"       => string(np.slopes),
        "slopelengths" => string(np.slopelengths),
        "values"       => string(np.values),
        "slopesbefore" => string(np.slopesbefore),
        "polystr"      => string(f)
    )

    println("Found new Newton polygon:")# $row")
    println("$np")
    push!(df, row)

    return df
end

function add_heights_patch!(c::SupabaseClient, heights::DataFrame. heights_table::String, atomic::Bool=false)
    server = fetch_table(c, heights_table)
    server_map = Dict(row.id => row.num for row in eachrow(server))

    for row in eachrow(heights)
        delta = row.num
        delta == 0 && continue

        current = get(server_map, row.id, 0)
        newval = current + delta

        url = "$(c.url)/rest/v1/$(heights_table)?id=eq.$(row.id)"
        payload = JSON.json(Dict("num" => newval))

        resp = HTTP.patch(
            url,
            headers = headers(c),
            body = payload
        )

        resp.status >= 300 && error(
            "PATCH failed for id=$(row.id): $(String(resp.body))"
        )
    end

    heights.num .= 0

    return heights
end

function add_heights_patch!(c::SupabaseClient, heights::DataFrame, heights_table::String, atomic::Bool=true)
    rows = [
        Dict("id" => row.id, "num" => row.num)
        for row in eachrow(heights)
        if row.num != 0
    ]

    isempty(rows) && return nothing

    url = "$(c.url)/rest/v1/rpc/add_heights_generic"
    payload = JSON.json(Dict(
        "target_table" => heights_table,
        "rows" => rows
    ))

    resp = HTTP.post(
        url,
        headers = [
            "apikey"        => c.key,
            "Authorization" => "Bearer $(c.key)",
            "Content-Type"  => "application/json"
        ],
        body = payload
    )

    if resp.status >= 300
        error("Supabase RPC failed: $(String(resp.body))")
    end

    heights.num .= 0
    return heights
end

############################################################
# PIPELINE
############################################################

function my_example(p)
    # R, (x1,x2,x3,x4) = GF(p)[:x1,:x2,:x3,:x4]
    R, (x1,x2,x3) = GF(p)[:x1,:x2,:x3]
    # f = x1^4 + x2^4 + x3^4 + x4^4 + 2x1*x2*x3*x4
    # penc(a) = x^6 + y^6 + z^6 + a*(x^5*y + y^5*x + x^5*z + z^5*x + y^5*z + z^5*y)

    # x^6 + 3*x^3*y^3 + y^6 + 3*y^3*z^3 + 4*z^6
    # f = x1^5 + x2^5 + x3^5 + x4^5
    # f

    # 6*x1^6 + 7*x1^5*x2 + 7*x1^5*x3 + x1^4*x2^2 + 4*x1^4*x2*x3 + x1^4*x3^2 + 10*x1^3*x2^3 + 4*x1^3*x2^2*x3 + 4*x1^3*x2*x3^2 + 10*x1^3*x3^3 + x1^2*x2^4 + 4*x1^2*x2^3*x3 + x1^2*x2^2*x3^2 + 4*x1^2*x2*x3^3 + x1^2*x3^4 + 7*x1*x2^5 + 4*x1*x2^4*x3 + 4*x1*x2^3*x3^2 + 4*x1*x2^2*x3^3 + 4*x1*x2*x3^4 + 7*x1*x3^5 + 6*x2^6 + 7*x2^5*x3 + x2^4*x3^2 + 10*x2^3*x3^3 + x2^2*x3^4 + 7*x2*x3^5 + 6*x3^6

    # 4*x^6 + 3*x^3*y^3 + y^6 + 3*y^3*z^3 + 4*z^6
    x1^8 + x2^8 + x3^8
end

function run_pipeline()

    start_time = time()

    heights = DataFrame(
        id  = collect(1:11),
        num = zeros(Int, 11)
    )

    println("Number of threads: $(Threads.nthreads())")

    println("Enter Supabase URL: ")
    url = "https://oigzppmsalxtinyuegpy.supabase.co"
    println("Enter Supabase API key (anon key): ")
    key = readline()

    println("Enter n: "); n = parse(Int, readline())
    println("Enter d: "); d = parse(Int, readline())
    println("Enter p: "); p = parse(Int, readline())
    println("Enter number of random samples N: "); N = parse(Int, readline())

    table = "n$(n)_d$(d)_np"
    println("Enter table name; leave blank for default $(table)"); t = readline()
    table = length(t) == 0 ? table : t

    heights_table = "heights"
    println("Enter heights table name; leave blank for default $(heights_table)"); t = readline()
    heights_table = length(t) == 0 ? heights_table : t

    client = SupabaseClient(url, key)

    println("\nFetching existing table…")
    df_old = fetch_table(client, table)
    println("Loaded $(nrow(df_old)) rows.\n")

    old_count = nrow(df_old)

    # RUN EXPERIMENT
    println("Running experiment…")
    # @time df = cpu_example_fast_random(n,d,p,N,df_old)
    # @time df = cpu_example_fast_example(n,d,p,N,my_example(p),df_old)
    # @time df = cpu_vector_fast_random(n,d,p,N,df_old)
    @time df = cpu_vector_fast_random_K3(n,d,p,N,df_old)
    # @time df = cpu_vector_fast_example(n,d,p,N,my_example(p),df_old)
    # @time df = cpu_vector_fast_symmetric(n,d,p,N,df_old)
    # @time df = cpu_vector_fast_weighted_example(n,d,p,N,my_example(p),df_old,num_monomials=2,nonzero_weights=true)

    # TO RUN WITH HEIGHTS
    # @time df, heights = func()

    # SAVE TO CSV
    filename = "$(table).csv"
    println("Saving full table to CSV: $filename")
    CSV.write(filename, df)

    # PUSH TO SUPABASE
    println("Pushing new rows to Supabase…")
    newrows = df[old_count+1:end, :]
    for row in eachrow(newrows)
        row_dict = Dict(names(newrows) .=> values(row))
        push_row(client, table, row_dict)
    end

    println("Pushing new heights to Supabase…")
    add_heights_patch!(client, heights, heights_table, atomic=false)

    println("Done.")

    println("Time taken: $(time() - start_time) seconds")
    
    return df[df.n .== n .&& df.d .== d .&& df.p .== p, :]
end

run_pipeline()
