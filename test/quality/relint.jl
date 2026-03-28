using ReLint
# to run do julia relint.jl FILEPATH

include("common.jl")

target = parse_target(ARGS)
files = target_julia_files(target)

println("ReLint loaded for TARGET=$(target)")
println("Julia files discovered: $(length(files))")

had_errors = false

for file in files
    println("\n==> $file")
    diags = ReLint.lint_file(file)
    for d in diags
        println(d)
    end
    had_errors |= !isempty(diags)
end

exit(had_errors ? 1 : 0)