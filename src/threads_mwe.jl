
struct MyStateStruct{T}
    A::T
    B::T
end

function do_work(s::MyStateStruct,index)
    println("$index, threadid=$(Threads.threadid()), A_id=$(objectid(s.A)), B_id=$(objectid(s.B))")
    if index == 1
        println("$index, threadid=$(Threads.threadid()), A_id=$(objectid(s.A)), B_id=$(objectid(s.B))")
    end
    C = s.A*s.B
    sum(C)
end

function new_random_struct(len)
    A = rand(1:100,len,len)
    B = rand(1:100,len,len)
    
    MyStateStruct(A,B)
end

function reset!(s::MyStateStruct)
    fill!(A,zero(eltype(A)))
    fill!(B,zero(eltype(A)))
end

function main()

    # create 
    nT = Threads.nthreads()
    state_channel = Channel{MyStateStruct{Matrix{Int}}}(nT)
    for i in 1:nT
        s = new_random_struct(500)
        println("$i, threadid=$(Threads.threadid()), A_id=$(objectid(s.A)), B_id=$(objectid(s.B))")
        put!(state_channel, s)
    end

    println("---")
    result = zeros(Int,30)

    Threads.@threads for i in 1:30

        s = take!(state_channel)
        println("$i, threadid=$(Threads.threadid()), A_id=$(objectid(s.A)), B_id=$(objectid(s.B))")

        do_work(s, i)
        put!(state_channel, s)

    end
    
end
