module Algebra

using StaticArrays
using StaticArrays: cross, dot
using ..Types

export evalNormal,
    evalNormalTriangle,
    evalNormalQuad,
    computeNormal,
    interpolateNormal,
    transformNormal,
    interpolateTriangle,
    interpolateQuad,
    norm,
    length,
    transformPoint,
    transformVector,
    transformDirection,
    unitVector,
    linInterp,
    sampleDisk

function evalNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if intersection.isTriangle
        return evalNormalTriangle(shape, intersection, frame)
    else
        return evalNormalQuad(shape, intersection, frame)
    end
end

function evalNormalTriangle(
    shape::Shape,
    intersection::Intersection,
    frame::Frame,
)
    (indexA, indexB, indexC) =
        @view shape.triangles[intersection.elementIndex, :]
    if (!isempty(shape.normals))
        normalA::SVec3f = SVec3f(@view shape.normals[indexA, :])
        normalB::SVec3f = SVec3f(@view shape.normals[indexB, :])
        normalC::SVec3f = SVec3f(@view shape.normals[indexC, :])
        normal = interpolateNormal(
            normalA,
            normalB,
            normalC,
            intersection.u,
            intersection.v,
        )
    else
        normal = computeNormal(
            SVec3f(@view shape.positions[indexA, :]),
            SVec3f(@view shape.positions[indexB, :]),
            SVec3f(@view shape.positions[indexC, :]),
        )
    end
    return transformNormal(frame, normal)
end

function evalNormalQuad(shape::Shape, intersection::Intersection, frame::Frame)
    (indexA, indexB, indexC, indexD) =
        @view shape.quads[intersection.elementIndex, :]
    if (!isempty(shape.normals))
        normalA::SVec3f = SVec3f(@view shape.normals[indexA, :])
        normalB::SVec3f = SVec3f(@view shape.normals[indexB, :])
        normalC::SVec3f = SVec3f(@view shape.normals[indexC, :])
        normalD::SVec3f = SVec3f(@view shape.normals[indexD, :])
        normal = interpolateNormal(
            normalA,
            normalB,
            normalC,
            normalD,
            intersection.u,
            intersection.v,
        )
    else
        normal = computeNormal(
            SVec3f(@view shape.positions[indexA, :]),
            SVec3f(@view shape.positions[indexB, :]),
            SVec3f(@view shape.positions[indexC, :]),
            SVec3f(@view shape.positions[indexD, :]),
        )
    end
    return transformNormal(frame, normal)
end

# computes the normal of a triangle
function computeNormal(pointA::SVec3f, pointB::SVec3f, pointC::SVec3f)
    return norm(cross(pointB .- pointA, pointC .- pointA))
end

# computes the normal of a quad 
function computeNormal(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
    pointD::SVec3f,
)
    return norm(
        computeNormal(pointA, pointB, pointD) +
        computeNormal(pointC, pointD, pointB),
    )
end

# interpolates the normals of a triangle
function interpolateNormal(
    normalA::SVec3f,
    normalB::SVec3f,
    normalC::SVec3f,
    u::Float32,
    v::Float32,
)
    return norm(interpolateTriangle(normalA, normalB, normalC, u, v))
end

# interpolates the normals of a quad
function interpolateNormal(
    normalA::SVec3f,
    normalB::SVec3f,
    normalC::SVec3f,
    normalD::SVec3f,
    u::Float32,
    v::Float32,
)
    return norm(interpolateQuad(normalA, normalB, normalC, normalD, u, v))
end
@inline function transformNormal(frame::Frame, v::SVec3f)
    return norm(transformVector(frame, v))
end

@inline function interpolateTriangle(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    u::Float32,
    v::Float32,
)::SVec3f
    return p0 * (1 - u - v) .+ p1 * u + p2 * v
end

@inline function interpolateQuad(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    p3::SVec3f,
    u::Float32,
    v::Float32,
)
    if (u + v <= 1)
        return interpolateTriangle(p0, p1, p3, u, v)
    else
        return interpolateTriangle(p2, p3, p1, 1 - u, 1 - v)
    end
end

@inline function norm(v::SVec3f)::SVec3f
    l = length(v)
    return (l != 0) ? v / l : v
end

@inline function length(v::SVec3f)::Float32
    return sqrt(dot(v, v))
end

@inline function transformPoint(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3] + frame.o
end

@inline function transformVector(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3]
end

@inline function transformDirection(frame::Frame, v::SVec3f)::SVec3f
    return norm(transformVector(frame, v))
end

@inline function unitVector(v::SVec3f)::SVec3f
    return v / length(v)
end

@inline function linInterp(a::SVec3f, b::SVec3f, weight::Float32)
    return a * (1 - weight) + b * weight
end

@inline function sampleDisk(u::Float32, v::Float32)::SVec2f
    r = sqrt(v)
    phi = 2 * pi * u
    return cos(phi) * r, sin(phi) * r
end

# end module
end