module MaterialsFunctions
using ..World
using ..Types
using ..Algebra

using StaticArrays: dot, cross

export sameHemisphere,
    reflectivityToEta,
    fresnelDielectric,
    fresnelConductor,
    microfacetDistribution,
    preciseMicrofacetDistribution,
    microfacetShadowing,
    sampleHemisphere,
    sampleHemisphereCos,
    sampleHemisphereCosPower,
    sampleMicrofacet,
    basisFromz

@inline function sameHemisphere(
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::Bool
    return dot(normal, outgoing) * dot(normal, incoming) >= 0
end

# @inline function fresnelSchlick(
#     specular::SVec3f,
#     normal::SVec3f,
#     outgoing::SVec3f,
# )::SVec3f
#     if specular == zeroSV3f
#         return zeroSV3f
#     end
#     cosine = dot(normal, outgoing)
#     specular .+
#     (1.0f0 .- specular) .* clamp(1.0f0 .- abs(cosine), 0.0f0, 1.0f0)^5.0f0
# end

@inline function reflectivityToEta(color::SVec3f)::SVec3f
    reflectivity::SVec3f = clamp.(color, 0.0f0, 0.99f0)
    return (1.0f0 .+ sqrt.(reflectivity)) / (1.0f0 .- sqrt.(reflectivity))
end

@inline function fresnelDielectric(
    eta::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::Float32
    cosw = abs(dot(normal, outgoing))

    sin2 = muladd(-cosw, cosw, 1.0f0)
    #sin2 = 1.0f0 - cosw * cosw
    eta2 = eta * eta

    cos2t = 1.0f0 - sin2 / eta2
    if cos2t < 0
        return 1.0f0
    end

    t0 = sqrt(cos2t)
    t1 = eta * t0
    t2 = eta * cosw

    rs = (cosw - t1) / (cosw + t1)
    rp = (t0 - t2) / (t0 + t2)

    #return (rs * rs + rp * rp) / 2.0f0
    return muladd(rs, rs, rp * rp) / 2.0f0 # fma seems faster but it may depends on architecture
end

@inline function fresnelConductor(
    eta::SVec3f,
    etak::SVec3f,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    cosw::Float32 = dot(normal, outgoing)
    if cosw <= 0
        return zeroSV3f
    end
    cosw = clamp(cosw, -1.0f0, 1.0f0)
    cos2::Float32 = cosw * cosw
    sin2::Float32 = clamp(1 - cos2, 0.0f0, 1.0f0)
    eta2::SVec3f = eta * eta
    etak2::SVec3f = etak * etak

    t0::SVec3f = eta2 - etak2 .- sin2
    a2plusb2::SVec3f = sqrt.(muladd.(t0, t0, 4.0f0 * eta2 * etak2))
    #a2plusb2::SVec3f = sqrt.(t0 * t0 + 4.0f0 * eta2 * etak2)
    t1::SVec3f = a2plusb2 .+ cos2
    a::SVec3f = sqrt.((a2plusb2 + t0) / 2.0f0)
    t2::SVec3f = 2.0f0 * a * cosw
    rs::SVec3f = (t1 - t2) / (t1 + t2)

    t3::SVec3f = muladd.(a2plusb2, cos2, sin2 * sin2)
    # t3::SVec3f = cos2 * a2plusb2 .+ sin2 * sin2
    t4::SVec3f = t2 * sin2
    rp::SVec3f = rs * (t3 - t4) / (t3 + t4)

    return (rp + rs) / 2.0f0
end

@inline function microfacetDistribution(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    ggx::Bool = true,
)::Float32
    cosine = dot(normal, halfway)
    if cosine <= 0
        return 0.0f0
    end
    roughness2 = roughness * roughness
    cosine2 = cosine * cosine
    if ggx
        return roughness2 / (
            pi *
            muladd(cosine2, roughness2, 1 - cosine2) *
            muladd(cosine2, roughness2, 1 - cosine2)
            #(cosine2 * roughness2 + 1 - cosine2) *
            #(cosine2 * roughness2 + 1 - cosine2)
        )
    else
        return exp((cosine2 - 1) / (roughness2 * cosine2)) /
               (pi * roughness2 * cosine2 * cosine2)
    end
end

@inline function preciseMicrofacetDistribution(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
)::Float32
    cosine = dot(normal, halfway)
    nxh = cross(normal, halfway)

    a = cosine * roughness
    k = roughness / muladd(a, a, dot(nxh, nxh))
    #k = roughness / (dot(nxh, nxh) + a * a)
    return k * k / pi
end

@inline function microfacetShadowing(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
    ggx::Bool = true,
)::Float32
    return microfacetShadowing1(roughness, normal, halfway, outgoing, ggx) *
           microfacetShadowing1(roughness, normal, halfway, incoming, ggx)
end

@inline function microfacetShadowing1(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    direction::SVec3f,
    ggx::Bool,
)::Float32
    cosine = dot(normal, direction)
    cosineh = dot(halfway, direction)
    if cosine * cosineh <= 0
        return 0.0f0
    end
    roughness2 = roughness * roughness
    cosine2 = cosine * cosine
    if ggx
        return 2.0f0 * abs(cosine) / (
            abs(cosine) +
            sqrt(muladd(-roughness2, cosine2, cosine2) + roughness2)
        )
        #(abs(cosine) + sqrt(cosine2 - roughness2 * cosine2 + roughness2))
    else
        ci = abs(cosine) / (roughness * sqrt(1 - cosine2))
        return ci < 1.6f0 ?
               muladd(3.535f0, ci, 2.181f0 * ci * ci) /
               (1.0f0 + muladd(2.276f0, ci, 2.577f0 * ci * ci)) : 1.0f0
        #    (3.535f0 * ci + 2.181f0 * ci * ci) /
        #    (1.0f0 + 2.276f0 * ci + 2.577f0 * ci * ci) : 1.0f0
    end
end

@inline function sampleHemisphere(normal::SVec3f)::SVec3f
    ruvx = rand(Float32)
    z = rand(Float32)  # ruvy
    r = sqrt(clamp(muladd(-z, z, 1.0f0), 0.0f0, 1.0f0))
    #r = sqrt(clamp(1 - z * z, 0.0f0, 1.0f0))
    phi = 2 * pi * ruvx
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

@inline function sampleHemisphereCos(normal::SVec3f)::SVec3f
    ruvx = rand(Float32)
    z = sqrt(rand(Float32))  # ruvy
    r = sqrt(muladd(-z, z, 1.0f0))
    #r = sqrt(1.0f0 - z * z)
    phi = 2.0f0 * pi * ruvx
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

@inline function sampleHemisphereCosPower(
    exponent::Float32,
    normal::SVec3f,
)::SVec3f
    ruxy = rand(Float32)
    ruvy = rand(Float32)
    z = ruvy^(1.0f0 / (exponent + 1.0f0))
    r = sqrt(muladd(-z, z, 1.0f0))
    #r = sqrt(1.0f0 - z * z)
    phi = 2.0f0 * pi * ruxy
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

@inline function sampleMicrofacet(
    roughness::Float32,
    normal::SVec3f,
    ggx::Bool = true,
)::SVec3f
    ruvx = rand(Float32)
    ruvy = rand(Float32)
    phi = 2.0f0 * pi * ruvx
    theta = 0.0f0
    if ggx
        theta = atan(roughness * sqrt(ruvy / (1.0f0 - ruvy)))
    else
        roughness2 = roughness * roughness
        theta = atan(sqrt(-roughness2 * log(1.0f0 - ruvy)))
    end
    localHalfVector =
        SVec3f(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta))
    return transformDirection(basisFromz(normal), localHalfVector)
end

#end module
end