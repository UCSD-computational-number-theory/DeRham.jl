using HTTP
using JSON
using CSV
using DataFrames
using DeRham

############################################################
# SUPABASE CLIENT
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

############################################################
# GET TABLE AS DATAFRAME
############################################################

function fetch_table(c::SupabaseClient, table::String)
    url = "$(c.url)/rest/v1/$table?select=*"
    resp = HTTP.get(url, headers=headers(c))

    rows = JSON.parse(String(resp.body))
    df = DataFrame(rows)
    return df
end

############################################################
# UPSERT A SINGLE ROW
############################################################

function push_row(c::SupabaseClient, table::String, row::Dict)
    baseurl = "$(c.url)/rest/v1/$table"
    payload = JSON.json(row)

    HTTP.post(
        baseurl,
        headers = [
            "apikey"        => c.key,
            "Authorization" => "Bearer $(c.key)",
            "Content-Type"  => "application/json",
            "Prefer"        => "resolution=merge-duplicates"
        ],
        body = payload
    )
end

############################################################
# ENCODERS
############################################################

encode_key(np) = JSON.json(np)
encode_value(f) = string(f)

############################################################
# MAIN EXPERIMENT LOOP
############################################################

function cpu_example_fast_random(n,d,p,N,df::DataFrame)
    T = d < n ? collect(0:d-1) : collect(0:n-1)
    lock = ReentrantLock()

    newrows = Vector{Dict}()

    Threads.@threads for i = 1:N
        f = DeRham.random_hypersurface(n,d,p)
        np = DeRham.newton_polygon(f, S=T, fastevaluation=true, algorithm=:naive)
        np === false && continue

        slopes       = np.slopes
        slopelengths = np.slopelengths
        values       = np.values
        slopesbefore = np.slopesbefore

        @lock lock begin
            # Check if an identical row exists
            same = filter(row ->
                row.n == n &&
                row.d == d &&
                row.p == p &&
                row.values == string(np.values),
            df)



            if nrow(same) == 0
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

                push!(newrows, row)
                push!(df, row)
            end
        end
    end

    return df, newrows
end

############################################################
# PIPELINE
############################################################

function run_pipeline()
    println("Enter Supabase URL: ")
    url = "https://oigzppmsalxtinyuegpy.supabase.co"
    println("Enter Supabase API key (anon key): ")
    key = readline()

    println("Enter n: "); n = parse(Int, readline())
    println("Enter d: "); d = parse(Int, readline())
    println("Enter p: "); p = parse(Int, readline())
    println("Enter number of random samples N: "); N = parse(Int, readline())

    client = SupabaseClient(url, key)
    table  = "derham"

    println("\nFetching existing table…")
    df = fetch_table(client, table)
    println("Loaded $(nrow(df)) rows.\n")

    println("Running experiment…")
    df, newrows = cpu_example_fast_random(n,d,p,N,df)

    println("New rows generated: ", length(newrows))

    ########################################################
    # SAVE FULL TABLE TO CSV
    ########################################################
    filename = "$(table).csv"
    println("Saving full table to CSV: $filename")
    CSV.write(filename, df)

    ########################################################
    # PUSH ONLY NEW ROWS TO SUPABASE
    ########################################################
    println("Pushing new rows to Supabase…")
    for row in newrows
        push_row(client, table, row)
    end

    println("Done.")
    return df
end

run_pipeline()
