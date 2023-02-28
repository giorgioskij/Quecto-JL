using BenchmarkTools
using StaticArrays

const SVec3 = SVector{3,Float32}

struct Vec3
    x::Float32
    y::Float32
    z::Float32
    v::SVec3
end

add(a::Vec3, b::Vec3)::Vec3 = Vec3(a.x + b.x, a.y + b.y, a.z + b.z, a.v .+ b.v)
mul(a::Vec3, b::Vec3)::Vec3 = Vec3(a.x * b.x, a.y * b.y, a.z * b.z, a.v * b.v)

function test(v1, v2)
    v3 = Vec3(0, 0, 0, SVec3([0, 0, 0]))
    for i = 1:10^7
        v3 = add(mul(v1, v2), v3)
    end
    return v3
end

v1 = Vec3(rand(), rand(), rand(), SVec3(rand(Float32, 3)))
v2 = Vec3(rand(), rand(), rand(), SVec3(rand(Float32, 3)))

test(v1, v2)
@benchmark v3 = test(v1, v2)

# println(v3)