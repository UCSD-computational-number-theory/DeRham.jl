using JET
using DeRham

target = isempty(ARGS) ? "." : ARGS[1]
println("Running JET report_package for DeRham (TARGET=$(target))")
report_package(DeRham; target_modules=(DeRham,))