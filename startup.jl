# DSMLP-friendly Julia startup script
# Initializes a writable user environment based on system-installed env

using Pkg
using InteractiveUtils

# Julia version environment path
const JVERSION = "v$(VERSION.major).$(VERSION.minor)"

# Paths
const SYS_ENV = "/usr/local/share/julia/environments/$(JVERSION)"
const USER_ROOT = "/data"
const USER_ENV = "/data/environments/$(JVERSION)"

# Ensure user depot path exists
if !isdir(USER_ENV)
    mkpath("/data/environments")
    cp(SYS_ENV, USER_ENV; force=true, recursive=true)
end

# Make sure depot is writable
pushfirst!(DEPOT_PATH, USER_ROOT)

# Activate the per-user environment
try
    Pkg.activate(USER_ENV; io=devnull)
catch e
    @warn "Could not activate user environment" exception=e
end
