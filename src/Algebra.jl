module Algebra

using StaticArrays
using StaticArrays: cross, dot
using ..Types

export transformNormal,
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
    transformRay,
    srgbToRgb,
    rgbToSrgb

@inline function transformNormal(
    frame::Frame,
    v::SVec3f,
    nonRigid = false,
)::SVec3f
    if nonRigid
        error("not implemented")
    else
        return norm(transformVector(frame, v))
    end
end

@inline function interpolateTriangle(
    p0::SVec2f,
    p1::SVec2f,
    p2::SVec2f,
    u::Float32,
    v::Float32,
)::SVec2f
    return p0 * (1 - u - v) + p1 * u + p2 * v
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
    p0::SVec2f,
    p1::SVec2f,
    p2::SVec2f,
    p3::SVec2f,
    u::Float32,
    v::Float32,
)::SVec2f
    if u + v <= 1
        return interpolateTriangle(p0, p1, p3, u, v)
    else
        return interpolateTriangle(p2, p3, p1, 1 - u, 1 - v)
    end
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
    return ifelse(l != 0, v / l, v)
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

@inline function linInterp(a::SVec4f, b::SVec4f, weight::Float32)::SVec4f
    return a * (1 - weight) .+ b * weight
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
    return dot(a.x, StaticArrays.cross(a.y, a.z))
end

@inline function adjoint(a::Mat3f)::Mat3f
    return transposeMat(
        Mat3f(
            StaticArrays.cross(a.y, a.z),
            StaticArrays.cross(a.z, a.x),
            StaticArrays.cross(a.x, a.y),
        ),
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

@inline function srgbToRgb(srgb::SVec4f)::SVec4f
    SVec4f(srgbToRgb(srgb[1]), srgbToRgb(srgb[2]), srgbToRgb(srgb[3]), srgb[4])
    # return SVec4f(srgbToRgb.(srgb)..., srgb[4])
end

@inline function srgbToRgb(srgb::SVec3f)::SVec3f
    # SVec3f(srgbToRgb(srgb[1]), srgbToRgb(srgb[2]), srgbToRgb(srgb[3]))
    return srgbToRgb.(srgb)
end

@inline function rgbToSrgb(rgb::SVec4f)::SVec4f
    return SVec4f(rgbToSrgb(rgb.x), rgbToSrgb(rgb.y), rgbToSrgb(rgb.z), rgb[4])
end

@inline function rgbToSrgb(rgb::Float32)::Float32
    # return (rgb <= 0.0031308f0) ? 12.92f0 * rgb :
    #        (1 + 0.055f0) * (rgb^(1 / 2.4f0)) - 0.055f0
    return (rgb <= 0.0031308f0) ? (12.92f0 * rgb) :
           (1 + 0.055f0) * fastPow(rgb, (1 / 2.4f0)) - 0.055f0
end
@inline function srgbToRgb(srgb::Float32)::Float32
    # srgb <= 0.04045 ? (srgb / 12.92f0) :
    # ((srgb + 0.055f0) / (1.0f0 + 0.055f0))^2.4f0
    return (srgb <= 0.04045) ? (srgb / 12.92f0) :
           fastPow((srgb + 0.055f0) / (1.0f0 + 0.055f0), 2.4f0)
end

# wtf? yeah you read it right, this can be slightly faster in julia. I know.
@inline function fastPow(a::Float32, b::Float32)
    return exp(log(a) * b)
end

# end module
end