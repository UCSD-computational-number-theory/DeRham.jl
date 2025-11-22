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

    if isempty(arr)
        println("No existing row. Returning empty dict.")
        return Dict{String,String}()
    end

    row = arr[1]

    # Construct dictionary of key → polynomial (as before)
    np = Dict(
        "slopes"       => row["slopes"],
        "slopelengths" => row["slopelengths"],
        "values"       => row["values"],
        "slopesbefore" => row["slopesbefore"]
    )

    key = JSON.json(np)
    val = row["polystr"]

    return Dict(key => val)
end

# Patch row's data column
function update_row(c::SupabaseClient, table::String, n::Int, d::Int, p::Int, dict::Dict)
    baseurl = "$(c.url)/rest/v1/$table"

    valid_entries = filter(((k,_),) -> startswith(k, "{"), dict)
    if isempty(valid_entries)
        println("No valid entries to update. Skipping.")
        return
    end
    
    (k, polystr) = first(valid_entries)
    parsed = JSON.parse(k)

    payload = JSON.json(Dict(
        "n" => n,
        "d" => d,
        "p" => p,
        "slopes"       => parsed["slopes"],
        "slopelengths" => parsed["slopelengths"],
        "values"       => parsed["values"],
        "slopesbefore" => parsed["slopesbefore"],
        "polystr"      => polystr
    ))

    HTTP.post(baseurl,
        headers = [
            "apikey" => c.key,
            "Authorization" => "Bearer $(c.key)",
            "Content-Type" => "application/json",
            "Prefer" => "resolution=merge-duplicates"
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
                resultsdict[encode_key(np)] = encode_value(f)
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
    println("Existing data: $existing_data")

    println("Running computation…")
    local_raw = cpu_example_fast_random(n, d, p, N, existing_data)
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

# create_sysimage(
#     [:HTTP, :JSON, :BitIntegers, :CSV, :Combinatorics, :DataFrames, :LFUDACache, :LRUCache, :LinearAlgebra, :Memoize, :OhMyThreads, :Oscar, :PProf];
#     sysimage_path="SysImage/CpuSysImage.so",
#     include_transitive_dependencies=false,
# )