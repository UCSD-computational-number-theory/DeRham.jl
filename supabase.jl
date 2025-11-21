using HTTP
using JSON
using DeRham

########## SIMPLE SUPABASE FUNCTIONS ##########

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

function get_row(c::SupabaseClient, table::String, n::Int, d::Int, p::Int)
    url = "$(c.url)/rest/v1/$table?n=eq.$n&d=eq.$d&p=eq.$p"
    resp = HTTP.get(url, headers=headers(c))
    arr = JSON.parse(String(resp.body))

    if length(arr) == 0
        # Row does not exist, make it
        println("No existing row. Creating a new one…")
        init_payload = JSON.json(Dict(
            "n" => n,
            "d" => d,
            "p" => p,
            "data" => Dict()
        ))
        post_url = "$(c.url)/rest/v1/$table"
        HTTP.post(post_url, headers=headers(c), body=init_payload)
        return Dict()
    end

    # Return existing data dictionary
    return arr[1]["data"]
end

# Patch row's data column
function update_row(c::SupabaseClient, table::String, n::Int, d::Int, p::Int, newdata::Dict)
    url = "$(c.url)/rest/v1/$table?n=eq.$n&d=eq.$d&p=eq.$p"

    payload = JSON.json(Dict("data" => newdata))

    HTTP.patch(url,
        headers = [
            "apikey" => c.key,
            "Authorization" => "Bearer $(c.key)",
            "Content-Type" => "application/json",
            "Prefer" => "return=representation"
        ],
        body = payload
    )
end

########## NEWTON POLYGON STORAGE ##########

encode_key(np) = JSON.json(np)
encode_value(f) = string(f)

function convert_results(raw::Dict)
    out = Dict{String,String}()
    for (np, f) in raw
        out[encode_key(np)] = encode_value(f)
    end
    return out
end

########## MAIN COMPUTE FUNCTION ##########

function cpu_example_fast_random(n,d,p,N,resultsdict)

    T = d < n ? collect(0:d-1) : collect(0:n-1)
    l = ReentrantLock()

    Threads.@threads for i = 1:N
        f = DeRham.random_hypersurface(n,d,p)
        np = DeRham.newton_polygon(f, S=T, fastevaluation=true, algorithm=:naive)

        if np != false
            @lock l begin
                if !haskey(resultsdict, np)
                    resultsdict[np] = f
                end
            end
        end
    end

    return resultsdict
end

########## PIPELINE ##########

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
    table = "derham"

    println("Fetching existing row…")
    existing_data = get_row(client, table, n, d, p)

    println("Running computation…")
    local_raw = cpu_example_fast_random(n, d, p, N, Dict())
    new_entries = convert_results(local_raw)

    println("Merging results…")
    merged = copy(existing_data)
    for (k,v) in new_entries
        merged[k] = v
    end

    println("Updating Supabase…")
    update_row(client, table, n, d, p, merged)

    println("Done.")
end

run_pipeline()