FROM alexpan5/qfsheight

RUN apt-get update && \  
    apt-get install --yes --no-install-recommends \ 
    curl ca-certificates vim git && \ 
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data

RUN git clone https://github.com/UCSD-computational-number-theory/DeRham.jl && \
    git clone https://github.com/UCSD-computational-number-theory/GPUFiniteFieldMatrices.jl