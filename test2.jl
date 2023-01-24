using StaticArrays
using BenchmarkTools

const Vec3 = SVector{3,Float64}

function test(v1, v2)
    v3 = Vec3(rand(), rand(), rand())
    for i = 1:10^7
        v3 = (v1 .* v2) .+ v3
    end
    return v3
end

v1 = Vec3(rand(), rand(), rand())
v2 = Vec3(rand(), rand(), rand())

test(v1, v2)
@benchmark v3 = test(v1, v2)
