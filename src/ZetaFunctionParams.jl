struct ZetaFunctionParams
    verbose::Int
    givefrobmat::Bool
    algorithm::Symbol
    termorder::Symbol
    vars_reversed::Bool
    fastevaluation::Bool
    always_use_bigints::Bool
    use_gpu::Bool
    use_threads::Bool
end

default_params() = ZetaFunctionParams(0,false,:default,:invlex,false,true,false,false,false)
