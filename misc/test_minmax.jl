using StaticArrays
using BenchmarkTools

const SVec3f = SVector{3,Float32}
const SVec4f = SVector{4,Float32}

@inbounds function fastMinMax(a::Float32, b::Float32)::Tuple{Float32,Float32}
    ifelse(Base.FastMath.lt_fast(b, a), (b, a), (a, b))
end

@inbounds function fastMinMaxVect(a::SVec3f, b::SVec3f)::Tuple{SVec3f,SVec3f}
    x, y, z =
        fastMinMax(a[1], b[1]), fastMinMax(a[2], b[2]), fastMinMax(a[3], b[3])
    return SVec3f(x[1], y[1], z[1]), SVec3f(x[2], y[2], z[2])
end

@inbounds function fastMinimumComponent(a::SVec4f)::Float32
    fastMin(fastMin(a[1], a[2]), fastMin(a[3], a[4]))
end

@inbounds function fastMaximumComponent(a::SVec4f)::Float32
    fastMax(fastMax(a[1], a[2]), fastMax(a[3], a[4]))
end

function fastMin(a::Float32, b::Float32)::Float32
    # ifelse(Base.FastMath.lt_fast(b, a), b, a)
    Base.FastMath.min_fast(a, b)
end

function fastMax(a::Float32, b::Float32)::Float32
    # ifelse(Base.FastMath.lt_fast(b, a), a, b)
    Base.FastMath.max_fast(a, b)
end

@inbounds function fastMaxVect(a::SVec3f, b::SVec3f)::SVec3f
    return SVec3f(fastMax(a[1], b[1]), fastMax(a[2], b[2]), fastMax(a[3], b[3]))
end

fastMinMaxVect(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand())) # to compile
@benchmark v3 = fastMinMaxVect(
    SVec3f(rand(), rand(), rand()),
    SVec3f(rand(), rand(), rand()),
)
# results: 21.409 ns +- 1.073 ns
# 0 allocation

fastMax.(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
@benchmark v4 =
    fastMax.(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results 20.580 ns +- 1.094 ns
# 0 allocation

fastMin.(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
@benchmark v5 =
    fastMin.(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results: 20.747 ns +- 1.845 ns
# 0 allocation

@benchmark v5 =
    map(fastMin, SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results: 20.957 ns +- 4.024 ns
# 0 allocation

@benchmark v6 =
    fastMaxVect(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results: 20.698 ns +- 1.565 ns
# 0 allocation

@benchmark v7 =
    min.(SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results: 19.249 ns +- 1.072 ns
# 0 allocation

@benchmark v8 =
    map(min, SVec3f(rand(), rand(), rand()), SVec3f(rand(), rand(), rand()))
# results: 19.320 ns +- 2.534 ns
# 0 allocation