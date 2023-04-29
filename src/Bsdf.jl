module Bsdf
using ..World
using ..Types
using ..Algebra
using ..MaterialsFunctions

using StaticArrays: dot, cross

export sampleBSDF, evalBSDF, evalDelta, sampleDelta

#####################################
# Sample BSDF
#####################################

@inline function sampleBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if (material.roughness == 0)
        return zeroSV3f
    end
    if material.type == "matte"
        return sampleMatte(normal, outgoing)
    elseif material.type == "glossy"
        return sampleGlossy(material.ior, material.roughness, normal, outgoing)
    elseif material.type == "reflective"
        return sampleReflective(material.roughness, normal, outgoing)
    elseif material.type == "transparent"
        return sampleTransparent(
            material.ior,
            material.roughness,
            normal,
            outgoing,
        )
    elseif material.type == "refractive"
        return sampleRefractive(
            material.ior,
            material.roughness,
            normal,
            outgoing,
        )
    elseif material.type == "subsurface"
        return sampleRefractive(
            material.ior,
            material.roughness,
            normal,
            outgoing,
        )
    else
        error("Unknown material type")
    end
end

@inline function sampleMatte(normal::SVec3f, outgoing::SVec3f)::SVec3f
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    incoming = sampleHemisphereCos(upNormal)
    return incoming
end

@inline function sampleGlossy(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    normal = dot(normal, outgoing) <= 0 ? -normal : normal
    if (rand(Float32) < fresnelDielectric(ior, normal, outgoing))
        halfway = sampleMicrofacet(roughness, normal)
        incoming = reflect(outgoing, halfway)
        if (!sameHemisphere(normal, outgoing, incoming))
            return zeroSV3f
        end
        return incoming
    else
        # incoming = sampleHemisphereCosPower(
        #     2.0f0 / (material.roughness * material.roughness),
        #     normal,
        # )
        incoming = sampleHemisphereCos(normal)
        return incoming
    end
end

@inline function sampleReflective(
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    # exponent = 2.0f0 / material.roughness^2
    # halfway = sampleHemisphereCosPower(exponent, normal)
    # incoming = reflect(outgoing, halfway)
    normal = dot(normal, outgoing) <= 0 ? -normal : normal
    halfway = sampleMicrofacet(roughness, normal)
    incoming = reflect(outgoing, halfway)
    if (!sameHemisphere(normal, outgoing, incoming))
        return zeroSV3f
    end
    return incoming
end

@inline function sampleTransparent(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    halfway = sampleMicrofacet(roughness, upNormal)
    if (rand(Float32) < fresnelDielectric(ior, halfway, outgoing))
        incoming = reflect(outgoing, halfway)
        if (!sameHemisphere(upNormal, outgoing, incoming))
            return zeroSV3f
        end
        return incoming
    else
        reflected = reflect(outgoing, halfway)
        incoming = -reflect(reflected, upNormal)
        if (sameHemisphere(upNormal, outgoing, incoming))
            return zeroSV3f
        end
        return incoming
    end
end

@inline function sampleRefractive(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)
    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    halfway = sampleMicrofacet(roughness, upNormal)
    if (
        rand(Float32) <
        fresnelDielectric(entering ? ior : (1.0f0 / ior), halfway, outgoing)
    )
        incoming = reflect(outgoing, halfway)
        if (!sameHemisphere(upNormal, outgoing, incoming))
            return zeroSV3f
        end
        return incoming
    else
        incoming = refract(outgoing, halfway, (entering ? (1.0f0 / ior) : ior))
        if (sameHemisphere(upNormal, outgoing, incoming))
            return zeroSV3f
        end
        return incoming
    end
end

#####################################
# Eval BSDF
#####################################

@inline function evalBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if (material.roughness == 0) #(dot(normal, incoming) * dot(normal, outgoing) <= 0) ||
        return zeroSV3f
    end

    if material.type == "matte"
        return evalMatte(material.color, normal, incoming, outgoing)

    elseif material.type == "glossy"
        return evalGlossy(
            material.color,
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "reflective"
        return evalReflective(
            material.color,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )

    elseif material.type == "transparent"
        return evalTransparent(
            material.color,
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "refractive"
        return evalRefractive(
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "subsurface"
        return evalRefractive(
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    else
        error("Unknown material type")
    end
end

@inline function evalMatte(
    color::SVec3f,
    normal::SVec3f,
    incoming::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if (dot(normal, incoming) * dot(normal, outgoing) <= 0)
        return zeroSV3f
    end
    radiance = color / pi * abs(dot(normal, incoming))
    return radiance
end

@inline function evalGlossy(
    color::SVec3f,
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return zeroSV3f
    end
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal

    F1 = fresnelDielectric(ior, upNormal, outgoing)
    halfway = norm(incoming + outgoing)
    F = fresnelDielectric(ior, halfway, incoming)
    D = microfacetDistribution(roughness, upNormal, halfway)
    # D = preciseMicrofacetDistribution(roughness, normal, halfway)
    G = microfacetShadowing(roughness, upNormal, halfway, outgoing, incoming)
    radiance =
        muladd.(
            color,
            (1 - F1) / pi * abs(dot(upNormal, incoming)),
            SVec3f(1, 1, 1) * F * D * G /
            (4.0f0 * dot(upNormal, outgoing) * dot(upNormal, incoming)) *
            abs(dot(upNormal, incoming)),
        )
    # color * (1 - F1) / pi * #2 *
    # abs(dot(upNormal, incoming)) +
    # SVec3f(1, 1, 1) * F * D * G /
    # (4.0f0 * dot(upNormal, outgoing) * dot(upNormal, incoming)) *
    # abs(dot(upNormal, incoming)) #^1.3
    return radiance
end

@inline function evalReflective(
    color::SVec3f,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return zeroSV3f
    end
    normal = dot(normal, outgoing) <= 0 ? -normal : normal

    halfway = norm(incoming + outgoing)
    F = fresnelConductor(reflectivityToEta(color), zeroSV3f, halfway, incoming)
    D = microfacetDistribution(roughness, normal, halfway)
    # D = preciseMicrofacetDistribution(roughness, normal, halfway)
    G = microfacetShadowing(roughness, normal, halfway, outgoing, incoming)
    radiance =
    #color / pi * 
        F * D * G / (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) * #* pi
        abs(dot(normal, incoming))
    return radiance
end

@inline function evalTransparent(
    color::SVec3f,
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        halfway = norm(incoming + outgoing)
        F = fresnelDielectric(ior, halfway, outgoing)
        D = microfacetDistribution(roughness, normal, halfway)
        G = microfacetShadowing(roughness, normal, halfway, outgoing, incoming)

        return SVec3f(1, 1, 1) * F * D * G /
               (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming))
    else
        reflected = reflect(-incoming, normal)
        halfway = norm(reflected + outgoing)
        F = fresnelDielectric(ior, halfway, outgoing)
        D = microfacetDistribution(roughness, normal, halfway)
        G = microfacetShadowing(roughness, normal, halfway, outgoing, reflected)

        return color * (1.0f0 - F) * D * G /
               (4.0f0 * dot(normal, outgoing) * dot(normal, reflected)) *
               (abs(dot(normal, reflected)))
    end
end

@inline function evalRefractive(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor::Float32 = entering ? ior : (1.0f0 / ior)
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        halfway = norm(incoming + outgoing)
        F = fresnelDielectric(relIor, halfway, outgoing)
        D = microfacetDistribution(roughness, upNormal, halfway)
        G = microfacetShadowing(
            roughness,
            upNormal,
            halfway,
            outgoing,
            incoming,
        )

        return SVec3f(1, 1, 1) * F * D * G /
               abs(4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming))
    else
        halfway =
            -norm(muladd.(relIor, incoming, outgoing)) *
            (entering ? 1.0f0 : -1.0f0)
        #-norm(relIor * incoming + outgoing) * (entering ? 1.0f0 : -1.0f0)
        F = fresnelDielectric(relIor, halfway, outgoing)
        D = microfacetDistribution(roughness, upNormal, halfway)
        G = microfacetShadowing(
            roughness,
            upNormal,
            halfway,
            outgoing,
            incoming,
        )

        return SVec3f(1, 1, 1) *
               abs(
                   (dot(outgoing, halfway) * dot(incoming, halfway)) /
                   (dot(outgoing, normal) * dot(incoming, normal)),
               ) *
               (1 - F) *
               D *
               G /
               (muladd(
            relIor,
            dot(halfway, incoming),
            dot(halfway, outgoing),
        )
        #relIor * dot(halfway, incoming) + dot(halfway, outgoing) # here if happer a comma from formatter delete it!
        )^2.0f0 * abs(dot(normal, incoming))
    end
end

#####################################
# Sample Delta
#####################################

@inline function sampleDelta(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
)
    if material.roughness != 0
        return zeroSV3f
    end

    if material.type == "reflective"
        return sampleDeltaReflective(normal, outgoing)
    elseif material.type == "transparent"
        return sampleDeltaTransparent(material.ior, normal, outgoing)
    elseif material.type == "refractive"
        return sampleDeltaRefractive(material.ior, normal, outgoing)
    elseif material.type == "volume" # not implemented for now
        sampleDeltaPassthrough(outgoing)
    else
        error("unknown material type")
    end
end

@inline function sampleDeltaReflective(normal::SVec3f, outgoing::SVec3f)::SVec3f
    normal = dot(normal, outgoing) <= 0 ? -normal : normal
    incoming = reflect(outgoing, normal)
    return incoming
end

@inline function sampleDeltaTransparent(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    if (rand(Float32) < fresnelDielectric(ior, upNormal, outgoing))
        return reflect(outgoing, upNormal)
    else
        return -outgoing
    end
end

@inline function sampleDeltaRefractive(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if (abs(ior - 1.0f0) < 1.0f-3)
        return -outgoing
    end

    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor = entering ? ior : (1.0f0 / ior)
    if (rand(Float32) < fresnelDielectric(relIor, upNormal, outgoing))
        return reflect(outgoing, upNormal)
    else
        return refract(outgoing, upNormal, 1.0f0 / relIor)
    end
end

@inline function sampleDeltaPassthrough(outgoing::SVec3f)::SVec3f
    return -outgoing
end

#####################################
# Eval Delta
#####################################

@inline function evalDelta(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if material.roughness != 0
        return zeroSV3f
    end

    if material.type == "reflective"
        return evalDeltaReflective(material.color, normal, outgoing, incoming)
    elseif material.type == "transparent"
        return evalDeltaTransparent(
            material.color,
            material.ior,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "refractive"
        return evalDeltaRefractive(material.ior, normal, outgoing, incoming)
    elseif material.type == "volume"
        return evalPassthrough(normal, outgoing, incoming)
    else
        error("unknown material type")
    end
end

@inline function evalDeltaReflective(
    color::SVec3f,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return zeroSV3f
    end
    normal = dot(normal, outgoing) <= 0 ? -normal : normal
    radiance =
    #material.color / pi * 
        fresnelConductor(reflectivityToEta(color), zeroSV3f, normal, outgoing)# / pi
    return radiance
end

@inline function evalDeltaTransparent(
    color::SVec3f,
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return SVec3f(1, 1, 1) * fresnelDielectric(ior, upNormal, outgoing)
    else
        return color * (1.0f0 - fresnelDielectric(ior, upNormal, outgoing))
    end
end

@inline function evalDeltaRefractive(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if (abs(ior - 1.0f0) < 1.0f-3)
        return dot(normal, incoming) * dot(normal, outgoing) <= 0 ?
               SVec3f(1, 1, 1) : zeroSV3f
    end

    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor = entering ? ior : (1.0f0 / ior)
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return SVec3f(1, 1, 1) * fresnelDielectric(relIor, upNormal, outgoing)
    else
        return SVec3f(1, 1, 1) *
               (1.0f0 / (relIor * relIor)) *
               (1.0f0 - fresnelDielectric(relIor, upNormal, outgoing))
    end
end

@inline function evalPassthrough(
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) >= 0
        return zeroSV3f
    else
        return SVec3f(1, 1, 1)
    end
end

#end module
end