# CUDA base image (with NVCC + toolkit + runtime)
FROM nvidia/cuda:12.1.1-devel-ubuntu20.04

# Basic tools
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --yes --no-install-recommends \
        curl ca-certificates git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Julia
ARG JULIA_VERSION=1.11.2
RUN curl -fsSL https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-$JULIA_VERSION-linux-x86_64.tar.gz \
    | tar -C /usr/local -xz --strip-components=1

# Create application directory
WORKDIR /app

# Copy your sysimage
COPY SysImage/DeRhamSysImage.so /usr/local/lib/julia/

# Optionally copy your package source if needed
# COPY DeRham/ /app/DeRham/

# Add a default runscript to test “using DeRham”
COPY supabase.jl /app/supabase.jl

# Use the sysimage by default
ENV JULIA_SYSIMAGE=/usr/local/lib/julia/DeRhamSysImage.so

# Start julia by default using this sysimage
ENTRYPOINT ["julia", "--sysimage=/usr/local/lib/julia/DeRhamSysImage.so"]
CMD ["julia", "supabase.jl"]
