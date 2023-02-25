module Types
using StaticArrays
using Images

export SVec4f, SVec3f, SVec2f, SVec3i, Frame, Mat3f, Ray, Triangle, Quad

const SVec4f = SVector{4,Float32}
const SVec3f = SVector{3,Float32}
const SVec2f = SVector{2,Float32}
const SVec3i = SVector{3,Int}
# const SVec4b = SVector{4,UInt8}

struct Frame
    x::SVec3f
    y::SVec3f
    z::SVec3f
    o::SVec3f

    Frame(x, y, z, o) = new(x, y, z, o)
    Frame(x, y, z) = new(x, y, z, SVec3f(0, 0, 0))
end

struct Mat3f
    x::SVec3f
    y::SVec3f
    z::SVec3f

    Mat3f(x, y, z) = new(x, y, z)
    Mat3f() = new(SVec3f(1, 0, 0), SVec3f(0, 1, 0), SVec3f(0, 0, 1))
end

struct Ray
    origin::SVec3f
    direction::SVec3f
    tmin::Float32
    tmax::Float32

    Ray(origin, direction, tmin, tmax) = new(origin, direction, tmin, tmax)
    Ray(origin, direction) = new(origin, direction, 0.0001f0, typemax(Float32))
end

struct Triangle
    x::SVec3f
    y::SVec3f
    z::SVec3f

    Triangle(x, y, z) = new(x, y, z)
end

struct Quad
    a::SVec3f
    b::SVec3f
    c::SVec3f
    d::SVec3f

    Quad(a, b, c, d) = new(a, b, c, d)
end

# end of module 
end