############################################################
# Base image
############################################################
ARG IMAGE=nvidia/cuda:12.1.1-devel-ubuntu20.04
FROM $IMAGE

############################################################
# System packages
############################################################
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --yes --no-install-recommends \
        curl ca-certificates git vim gcc g++ make && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

############################################################
# Install Julia (system install)
############################################################
ARG JULIA_RELEASE=1.12
ARG JULIA_VERSION=1.12.2

RUN curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_RELEASE}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz | \
    tar -C /usr/local -x -z --strip-components=1 -f -

############################################################
# System-wide Julia depot (root install)
############################################################
ENV JULIA_DEPOT_PATH=/usr/local/share/julia

RUN mkdir -p /usr/local/share/julia/environments/v${JULIA_RELEASE} && \
    chmod -R 0777 /usr/local/share/julia

############################################################
# Copy DeRham + GPUFiniteFieldMatrices into system depot env
############################################################
COPY Project.toml Manifest.toml supabase.jl /usr/local/share/julia/environments/v${JULIA_RELEASE}/
COPY src /usr/local/share/julia/environments/v${JULIA_RELEASE}/src

# External dependency
COPY --from=gf / /usr/local/share/julia/environments/v${JULIA_RELEASE}/GPUFiniteFieldMatrices

############################################################
# Instantiate & precompile system-wide packages
############################################################
RUN julia -e 'using Pkg; \
    Pkg.activate("/usr/local/share/julia/environments/v'${JULIA_RELEASE}'"); \
    Pkg.develop(path="/usr/local/share/julia/environments/v'${JULIA_RELEASE}'/GPUFiniteFieldMatrices"); \
    Pkg.add("HTTP"); \
    Pkg.add("JSON"); \
    Pkg.add("DataFrames")'

############################################################
# Runtime User Depot
############################################################
RUN mkdir -m 0777 /data
ENV JULIA_HISTORY=/data/logs/repl_history.jl

# Startup script to initialize user depot from system depot
COPY startup.jl /usr/local/share/julia/config/

############################################################
# Final runtime depot path (user depot + system depot)
############################################################
ENV JULIA_DEPOT_PATH=/data:/usr/local/share/julia
