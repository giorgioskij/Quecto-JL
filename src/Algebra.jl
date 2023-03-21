module Algebra

using StaticArrays
using StaticArrays: cross, dot
using ..Types
import Base.*
import Base./
using LoopVectorization
using Images

export transformNormal,
    interpolateTriangle,
    interpolateQuad,
    interpolateLine,
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
    rgbToSrgb,
    reflect,
    refract,
    basisFromz,
    orthonormalize

@inline function *(a::SVector, b::SVector)::SVector
    # Why this is faster than the broadcast variant?
    # on my PC it's seems broadcast is almost the same speed
    @turbo map(*, a, b)
    # why here we can't use vmap instead of map? It outputs monocolored image
end

@inline function /(a::SVector, b::SVector)::SVector
    # Why this is faster than the broadcast variant?
    # on my PC it's seems broadcast is almost the same speed
    @turbo map(/, a, b)
    # why here we can't use vmap instead of map? It outputs monocolored image
end

@inline @fastmath function transformNormal(
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

# TODO: merge all this duplicates is a more generic type
@inline @fastmath function interpolateTriangle(
    p0::SVec2f,
    p1::SVec2f,
    p2::SVec2f,
    u::Float32,
    v::Float32,
)::SVec2f
    return muladd.(p0, (1 - u - v), muladd.(p1, u, p2 * v))
    #return p0 * (1 - u - v) + p1 * u + p2 * v
end

@inline @fastmath function interpolateTriangle(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    u::Float32,
    v::Float32,
)::SVec3f
    return muladd.(p0, (1 - u - v), muladd.(p1, u, p2 * v))
    #return p0 * (1 - u - v) .+ p1 * u + p2 * v
end

@inline @fastmath function interpolateQuad(
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

@inline @fastmath function interpolateQuad(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    p3::SVec3f,
    u::Float32,
    v::Float32,
)::SVec3f
    if (u + v <= 1)
        return interpolateTriangle(p0, p1, p3, u, v)
    else
        return interpolateTriangle(p2, p3, p1, 1 - u, 1 - v)
    end
end

@inline @fastmath function interpolateLine(p0::SVec3f, p1::SVec3f, u::Float32)
    return p0 * (1 - u) + p1 * u
    # return muladd.(p0, (1 - u), p1 * u)
end

@inline @fastmath function interpolateLine(p0::SVec2f, p1::SVec2f, u::Float32)
    return p0 * (1 - u) + p1 * u
    # return muladd.(p0, (1 - u), p1 * u)
end

@inline @fastmath function norm(v::SVec3f)::SVec3f
    l = length(v)
    return ifelse(l != 0, v / l, v)
end

@inline @fastmath function length(v::SVec3f)::Float32
    return sqrt(dot(v, v))
end

@inline @inbounds @fastmath function transformPoint(
    frame::Frame,
    v::SVec3f,
)::SVec3f
    return muladd.(
        frame.x,
        v[1],
        muladd.(frame.y, v[2], muladd.(frame.z, v[3], frame.o)),
    )
    #return frame.x * v[1] + frame.y * v[2] + frame.z * v[3] + frame.o
end

@inline @inbounds @fastmath function transformVector(
    frame::Frame,
    v::SVec3f,
)::SVec3f
    return muladd.(frame.x, v[1], muladd.(frame.y, v[2], frame.z * v[3]))
    # return frame.x * v[1] + frame.y * v[2] + frame.z * v[3]
end

@inline @inbounds @fastmath function transformVector(
    a::Mat3f,
    v::SVec3f,
)::SVec3f
    return muladd.(a.x, v[1], muladd.(a.y, v[2], a.z * v[3]))
    #return a.x * v[1] + a.y * v[2] + a.z * v[3]
end

@inline function transformDirection(frame::Frame, v::SVec3f)::SVec3f
    return norm(transformVector(frame, v))
end

@inline function transformDirection(a::Mat3f, v::SVec3f)::SVec3f
    return norm(transformVector(a, v))
end

@inline @fastmath function unitVector(v::SVec3f)::SVec3f
    return v / length(v)
end

@inline @fastmath function linInterp(
    a::SVec4f,
    b::SVec4f,
    weight::Float32,
)::SVec4f
    return muladd.(a, (1 - weight), b * weight)
    # return a * (1 - weight) .+ b * weight
end

@inline @fastmath function linInterp(
    a::SVec3f,
    b::SVec3f,
    weight::Float32,
)::SVec3f
    a = map(clamp01nan, a)
    b = map(clamp01nan, b)

    return muladd.(a, (1 - weight), b * weight)
    #return a * (1 - weight) + b * weight
end

@inline @fastmath function sampleDisk(u::Float32, v::Float32)::SVec2f
    r = sqrt(v)
    phi = 2 * pi * u
    return cos(phi) * r, sin(phi) * r
end

@inline @inbounds @fastmath function determinant(a::Mat3f)::Float32
    return dot(a.x, StaticArrays.cross(a.y, a.z))
end

@inline @inbounds @fastmath function adjoint(a::Mat3f)::Mat3f
    return transposeMat(
        Mat3f(
            StaticArrays.cross(a.y, a.z),
            StaticArrays.cross(a.z, a.x),
            StaticArrays.cross(a.x, a.y),
        ),
    )
end

@inline @fastmath function inverse(a::Mat3f)::Mat3f
    return matMulFloat(adjoint(a), (1 / determinant(a)))
end

@inline @fastmath function inverse(frame::Frame, nonRigid::Bool = false)::Frame
    if nonRigid
        minv = inverse(rotation(frame))
        return makeFrame(minv, -(matMulVec(minv, frame.o)))
    else
        minv::Mat3f = transposeMat(rotation(frame))
        return makeFrame(minv, -(matMulVec(minv, frame.o)))
    end
end

@inline @inbounds function rotation(frame::Frame)::Mat3f
    return Mat3f(frame.x, frame.y, frame.z)
end

@inline @inbounds function transposeMat(mat::Mat3f)::Mat3f
    return Mat3f(
        SVec3f(mat.x[1], mat.y[1], mat.z[1]),
        SVec3f(mat.x[2], mat.y[2], mat.z[2]),
        SVec3f(mat.x[3], mat.y[3], mat.z[3]),
    )
end

@inline @inbounds function makeFrame(m::Mat3f, t::SVec3f)::Frame
    return Frame(m.x, m.y, m.z, t)
end

@inline @inbounds @fastmath function matMulVec(a::Mat3f, b::SVec3f)::SVec3f
    return SVec3f(muladd.(a.x, b[1], muladd.(a.y, b[2], a.z * b[3])))
    #return SVec3f(a.x * b[1] + a.y * b[2] + a.z * b[3])
end

@inline @inbounds @fastmath function matMulFloat(a::Mat3f, b::Float32)::Mat3f
    return Mat3f(a.x * b, a.y * b, a.z * b)
end

@inline @inbounds function xyz(a::SVec4f)::SVec3f
    return SVec3f(a[1], a[2], a[3])
end

@inline @fastmath function transformRay(frame::Frame, ray::Ray)::Ray
    return Ray(
        transformPoint(frame, ray.origin),
        transformVector(frame, ray.direction),
        ray.tmin,
        ray.tmax,
    )
end

@inline @inbounds function srgbToRgb(srgb::SVec4f)::SVec4f
    SVec4f(srgbToRgb(srgb[1]), srgbToRgb(srgb[2]), srgbToRgb(srgb[3]), srgb[4])
    # return SVec4f(srgbToRgb.(srgb)..., srgb[4])
end

@inline function srgbToRgb(srgb::SVec3f)::SVec3f
    # SVec3f(srgbToRgb(srgb[1]), srgbToRgb(srgb[2]), srgbToRgb(srgb[3]))
    return srgbToRgb.(srgb)
end

@inline @inbounds function rgbToSrgb(rgb::SVec4f)::SVec4f
    return SVec4f(rgbToSrgb(rgb.x), rgbToSrgb(rgb.y), rgbToSrgb(rgb.z), rgb[4])
end

@inline @fastmath function rgbToSrgb(rgb::Float32)::Float32
    # return (rgb <= 0.0031308f0) ? 12.92f0 * rgb :
    #        (1 + 0.055f0) * (rgb^(1 / 2.4f0)) - 0.055f0
    return ifelse(
        rgb <= 0.0031308f0,
        12.92f0 * rgb,
        muladd(1.055f0, fastPow(rgb, (1 / 2.4f0)), -0.055f0),
        #(1 + 0.055f0) * fastPow(rgb, (1 / 2.4f0)) - 0.055f0,
    )
end

@inline @fastmath function srgbToRgb(srgb::Float32)::Float32
    # srgb <= 0.04045 ? (srgb / 12.92f0) :
    # ((srgb + 0.055f0) / (1.0f0 + 0.055f0))^2.4f0
    return ifelse(
        srgb <= 0.04045f0,
        srgb / 12.92f0,
        fastPow((srgb + 0.055f0) / 1.055f0, 2.4f0),
        #fastPow((srgb + 0.055f0) / (1.0f0 + 0.055f0), 2.4f0),
    )
end

# wtf? yeah you read it right, this can be slightly faster in julia. I know.
@inline @fastmath function fastPow(a::Float32, b::Float32)::Float32
    return exp(log(a) * b)
end

@inline @fastmath function reflect(w::SVec3f, n::SVec3f)::SVec3f
    return muladd.(2.0f0 * dot(n, w), n, -w)
    #return -w + 2.0f0 * dot(n, w) * n
end

@inline @fastmath function refract(
    w::SVec3f,
    n::SVec3f,
    invEta::Float32,
)::SVec3f
    cosine = dot(n, w)
    k = muladd(invEta * invEta, muladd(cosine, cosine, -1.0f0), 1.0f0)
    #k = 1.0f0 + invEta * invEta * (cosine * cosine - 1.0f0)
    if k < 0
        return zeroSV3f
    end
    return muladd.(-w, invEta, muladd(invEta, cosine, -sqrt(k)) * n)
    #return -w * invEta + (invEta * cosine - sqrt(k)) * n
end

# Constructs a basis from a direction
@inline function basisFromz(v::SVec3f)
    z = norm(v)
    sign = copysign(1.0f0, z.z)
    a = -1.0f0 / (sign + z.z)
    b = z.x * z.y * a
    x = SVec3f(muladd(sign * z.x * z.x, a, 1.0f0), sign * b, -sign * z.x)
    #x = SVec3f(1.0f0 + sign * z.x * z.x * a, sign * b, -sign * z.x)
    y = SVec3f(b, muladd(z.y * z.y, a, sign), -z.y)
    #y = SVec3f(b, sign + z.y * z.y * a, -z.y)
    return Mat3f(x, y, z)
end

@inline function orthonormalize(a::SVec3f, b::SVec3f)::SVec3f
    return norm(a - b * dot(a, b))
end

@inline function lineTangent(p0::SVec3f, p1::SVec3f)::SVec3f
    return normalize(p1 - p0)
end

# end module
end