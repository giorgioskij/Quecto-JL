using BenchmarkTools

struct CustomType
    a::Float32
    b::Int
end

function composeList()
    n = 10_000_000
    v = Vector{CustomType}(undef, 0)
    v = Vector{CustomType}(undef, n)
    for i = 1:n
        @fastmath @inbounds v[i] = (CustomType(i / 5, i))
    end
    return v
end

function composeList2()
    v = Vector{CustomType}(undef, 0)
    n = 10_000_000
    sizehint!(v, n)
    for i = 1:n
        @fastmath @inbounds push!(v, CustomType(i / 5, i))
    end
    return v
end