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
    sampleDisk,
    inverse,
    rotation,
    transposeMat,
    makeFrame,
    matMulVec,
    xyz,
    transformRay

function evalNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if isempty(shape.normals)
        return computeNormal(shape, intersection, frame)
    end

    if !isempty(shape.triangles)
        (indexA, indexB, indexC) =
            @view shape.triangles[intersection.elementIndex, :]

        normalA = SVec3f(@view shape.normals[indexA, :])
        normalB = SVec3f(@view shape.normals[indexB, :])
        normalC = SVec3f(@view shape.normals[indexC, :])
        return transformNormal(
            frame,
            norm(
                interpolateTriangle(
                    normalA,
                    normalB,
                    normalC,
                    intersection.u,
                    intersection.v,
                ),
            ),
        )

    elseif !isempty(shape.quads)
        indexA, indexB, indexC, indexD =
            @view shape.quads[intersection.elementIndex, :]

        normalA = SVec3f(@view shape.normals[indexA, :])
        normalB = SVec3f(@view shape.normals[indexB, :])
        normalC = SVec3f(@view shape.normals[indexC, :])
        normalD = SVec3f(@view shape.normals[indexD, :])
        return transformNormal(
            frame,
            norm(
                interpolateQuad(
                    normalA,
                    normalB,
                    normalC,
                    normalD,
                    intersection.u,
                    intersection.v,
                ),
            ),
        )
    else
        error("Only triangles and quads right now")
    end
end

function computeNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if !isempty(shape.triangles)
        (indexA, indexB, indexC) =
            @view shape.triangles[intersection.elementIndex, :]
        pointA = SVec3f(@view shape.positions[indexA, :])
        pointB = SVec3f(@view shape.positions[indexB, :])
        pointC = SVec3f(@view shape.positions[indexC, :])
        return transformNormal(
            frame,
            computeTriangleNormal(pointA, pointB, pointC),
        )
    elseif !isempty(shape.quads)
        indexA, indexB, indexC, indexD =
            @view shape.quads[intersection.elementIndex, :]
        pointA = SVec3f(@view shape.positions[indexA, :])
        pointB = SVec3f(@view shape.positions[indexB, :])
        pointC = SVec3f(@view shape.positions[indexC, :])
        pointD = SVec3f(@view shape.positions[indexD, :])

        return transformNormal(
            frame,
            computeQuadNormal(pointA, pointB, pointC, pointD),
        )
    else
        error("Only triangles and quads right now")
    end
end

# function evalNormalTriangle(
#     shape::Shape,
#     intersection::Intersection,
#     frame::Frame,
# )
#     (indexA, indexB, indexC) =
#         @view shape.triangles[intersection.elementIndex, :]
#     if (!isempty(shape.normals))
#         normalA::SVec3f = SVec3f(@view shape.normals[indexA, :])
#         normalB::SVec3f = SVec3f(@view shape.normals[indexB, :])
#         normalC::SVec3f = SVec3f(@view shape.normals[indexC, :])
#         normal = interpolateNormal(
#             normalA,
#             normalB,
#             normalC,
#             intersection.u,
#             intersection.v,
#         )
#     else
#         normal = computeTriangleNormal(
#             SVec3f(@view shape.positions[indexA, :]),
#             SVec3f(@view shape.positions[indexB, :]),
#             SVec3f(@view shape.positions[indexC, :]),
#         )
#     end
#     return transformNormal(frame, normal)
# end

# function evalNormalQuad(shape::Shape, intersection::Intersection, frame::Frame)
#     (indexA, indexB, indexC, indexD) =
#         @view shape.quads[intersection.elementIndex, :]
#     if (!isempty(shape.normals))
#         normalA::SVec3f = SVec3f(@view shape.normals[indexA, :])
#         normalB::SVec3f = SVec3f(@view shape.normals[indexB, :])
#         normalC::SVec3f = SVec3f(@view shape.normals[indexC, :])
#         normalD::SVec3f = SVec3f(@view shape.normals[indexD, :])
#         normal = interpolateNormal(
#             normalA,
#             normalB,
#             normalC,
#             normalD,
#             intersection.u,
#             intersection.v,
#         )
#     else
#         normal = computeQuadNormal(
#             SVec3f(@view shape.positions[indexA, :]),
#             SVec3f(@view shape.positions[indexB, :]),
#             SVec3f(@view shape.positions[indexC, :]),
#             SVec3f(@view shape.positions[indexD, :]),
#         )
#     end
#     return transformNormal(frame, normal)
# end

# computes the normal of a triangle
@inline function computeTriangleNormal(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
)
    return norm(cross(pointB .- pointA, pointC .- pointA))
end

# computes the normal of a quad 
function computeQuadNormal(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
    pointD::SVec3f,
)
    return norm(
        computeTriangleNormal(pointA, pointB, pointD) +
        computeTriangleNormal(pointC, pointD, pointB),
    )
end

# # interpolates the normals of a triangle
# function interpolateNormal(
#     normalA::SVec3f,
#     normalB::SVec3f,
#     normalC::SVec3f,
#     u::Float32,
#     v::Float32,
# )
#     return norm(interpolateTriangle(normalA, normalB, normalC, u, v))
# end

# # interpolates the normals of a quad
# function interpolateNormal(
#     normalA::SVec3f,
#     normalB::SVec3f,
#     normalC::SVec3f,
#     normalD::SVec3f,
#     u::Float32,
#     v::Float32,
# )
#     return norm(interpolateQuad(normalA, normalB, normalC, normalD, u, v))
# end

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
    if isnan(a.x)
        a = SVec3f(0, a.y, a.z)
    end
    if isnan(a.y)
        a = SVec3f(a.x, 0, a.z)
    end
    if isnan(a.z)
        a = SVec3f(a.x, a.y, 0)
    end
    if isnan(b.x)
        b = SVec3f(0, b.y, b.z)
    end
    if isnan(b.y)
        b = SVec3f(b.x, 0, b.z)
    end
    if isnan(b.z)
        b = SVec3f(b.x, b.y, 0)
    end

    return a * (1 - weight) + b * weight
end

@inline function sampleDisk(u::Float32, v::Float32)::SVec2f
    r = sqrt(v)
    phi = 2 * pi * u
    return cos(phi) * r, sin(phi) * r
end

@inline function determinant(a::Mat3f)::Float32
    return dot(a.x, cross(a.y, a.z))
end

@inline function adjoint(a::Mat3f)::Mat3f
    return transposeMat(
        Mat3f(cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)),
    )
end

@inline function inverse(a::Mat3f)::Mat3f
    return matMulFloat(adjoint(a), (1 / determinant(a)))
end

@inline function inverse(frame::Frame, nonRigid::Bool = false)::Frame
    if nonRigid
        minv = inverse(rotation(frame))
        return makeFrame(minv, -(matMulVec(minv, frame.o)))
    else
        minv::Mat3f = transposeMat(rotation(frame))
        return makeFrame(minv, -(matMulVec(minv, frame.o)))
    end
end

@inline function rotation(frame::Frame)::Mat3f
    return Mat3f(frame.x, frame.y, frame.z)
end

@inline function transposeMat(mat::Mat3f)::Mat3f
    return Mat3f(
        SVec3f(mat.x[1], mat.y[1], mat.z[1]),
        SVec3f(mat.x[2], mat.y[2], mat.z[2]),
        SVec3f(mat.x[3], mat.y[3], mat.z[3]),
    )
end

@inline function makeFrame(m::Mat3f, t::SVec3f)::Frame
    return Frame(m.x, m.y, m.z, t)
end

@inline function matMulVec(a::Mat3f, b::SVec3f)::SVec3f
    return SVec3f(a.x * b[1] + a.y * b[2] + a.z * b[3])
end

@inline function matMulFloat(a::Mat3f, b::Float32)::Mat3f
    return Mat3f(a.x * b, a.y * b, a.z * b)
end

@inline function xyz(a::SVec4f)::SVec3f
    return SVec3f(a[1], a[2], a[3])
end

@inline function transformRay(frame::Frame, ray::Ray)::Ray
    return Ray(
        transformPoint(frame, ray.origin),
        transformVector(frame, ray.direction),
        ray.tmin,
        ray.tmax,
    )
end

# end module
end