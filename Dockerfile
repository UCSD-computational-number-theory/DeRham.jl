ARG IMAGE=nvidia/cuda:12.1.1-devel-ubuntu20.04
FROM $IMAGE

# ---------------------------------------------------------
# System tools
# ---------------------------------------------------------
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        curl ca-certificates git vim gcc g++ make && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------
# Install Julia (Linux build)
# ---------------------------------------------------------
ARG JULIA_RELEASE=1.12
ARG JULIA_VERSION=1.12.2

RUN curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_RELEASE}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz \
    | tar -C /usr/local -x -z --strip-components=1 -f -

# ---------------------------------------------------------
# Prepare depot
# ---------------------------------------------------------
ENV JULIA_DEPOT_PATH=/usr/local/share/julia

# ---------------------------------------------------------
# Copy DeRham project
# ---------------------------------------------------------
WORKDIR /opt/DeRham
COPY Project.toml .
COPY Manifest.toml .
COPY src/ src/

# ---------------------------------------------------------
# Bring GPUFiniteFieldMatrices from external context
# ---------------------------------------------------------
COPY --from=gf / /opt/GPUFiniteFieldMatrices/

# ---------------------------------------------------------
# Instantiate full environment (Linux artifacts!)
# ---------------------------------------------------------
RUN julia -e 'using Pkg; \
    Pkg.activate("/opt/DeRham"); \
    # Overwrite any existing path dependency from host
    Pkg.develop(path="/opt/GPUFiniteFieldMatrices"); \
    Pkg.resolve(); \
    Pkg.instantiate(); \
    Pkg.add("HTTP"); \
    Pkg.add("JSON");'

# ---------------------------------------------------------
# Runtime settings
# ---------------------------------------------------------

WORKDIR /opt/DeRham
