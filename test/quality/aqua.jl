using Aqua
using DeRham

target = isempty(ARGS) ? "." : ARGS[1]
println("Running Aqua for module DeRham (TARGET=$(target))")
Aqua.test_all(DeRham; stale_deps=false, deps_compat=false)