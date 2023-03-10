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
        return SVec3f(0, 0, 0)
    end
    if material.type == "matte"
        return sampleMatte(normal)
    elseif material.type == "glossy"
        return sampleGlossy(material.ior, material.roughness, normal, outgoing)
    elseif material.type == "reflective"
        return sampleReflective(material.roughness, normal, outgoing)
        # TODO: implement other materials
        # transparent
        # refractive
        # subsurface
    else
        error("Unknown material type")
    end
end

@inline function sampleMatte(normal::SVec3f)::SVec3f
    incoming = sampleHemisphereCos(normal)
    return incoming
end

@inline function sampleGlossy(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if (rand(Float32) < fresnelDielectric(ior, normal, outgoing))
        halfway = sampleMicrofacet(roughness, normal)
        incoming = reflect(outgoing, halfway)
        if (!sameHemisphere(normal, outgoing, incoming))
            return SVec3f(0, 0, 0)
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
    halfway = sampleMicrofacet(roughness, normal)
    incoming = reflect(outgoing, halfway)
    if (!sameHemisphere(normal, outgoing, incoming))
        return SVec3f(0, 0, 0)
    end
    return incoming
end

#####################################
# Eval BSDF
#####################################

function evalBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if (dot(normal, incoming) * dot(normal, outgoing) <= 0) ||
       (material.roughness == 0)
        return SVec3f(0, 0, 0)
    end

    if material.type == "matte"
        return evalMatte(material.color, normal, incoming)

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

        #TODO: implement other materials
        # transparent
        # refractive
        # subsurface
    else
        error("Unknown material type")
    end
end

@inline function evalMatte(
    color::SVec3f,
    normal::SVec3f,
    incoming::SVec3f,
)::SVec3f
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
        return SVec3f(0, 0, 0)
    end

    F1 = fresnelDielectric(ior, normal, outgoing)
    halfway = norm(incoming + outgoing)
    F = fresnelDielectric(ior, halfway, incoming)
    D = microfacetDistribution(roughness, normal, halfway)
    # D = preciseMicrofacetDistribution(roughness, normal, halfway)
    G = microfacetShadowing(roughness, normal, halfway, outgoing, incoming)
    radiance =
        color * (1 - F1) / pi * #2 *
        abs(dot(normal, incoming)) +
        SVec3f(1, 1, 1) * F * D * G /
        (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
        abs(dot(normal, incoming)) #^1.3
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
        return SVec3f(0, 0, 0)
    end

    halfway = norm(incoming + outgoing)
    F = fresnelConductor(
        reflectivityToEta(color),
        SVec3f(0, 0, 0),
        halfway,
        incoming,
    )
    D = microfacetDistribution(roughness, normal, halfway)
    # D = preciseMicrofacetDistribution(roughness, normal, halfway)
    G = microfacetShadowing(roughness, normal, halfway, outgoing, incoming)
    radiance =
    #color / pi * 
        F * D * G / (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) * #* pi
        abs(dot(normal, incoming))
    return radiance
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
        return SVec3f(0, 0, 0)
    end

    if material.type == "reflective"
        return sampleDeltaReflective(normal, outgoing)
    elseif material.type == "transparent"
        return sampleDeltaTransparent(material.ior, normal, outgoing)
    elseif material.type == "refractive"
        return sampleDeltaRefractive(material.ior, normal, outgoing)
        # elseif material.type == "volumetric" # not implemented for now
    else
        error("unknown material type")
    end
end

@inline function sampleDeltaReflective(normal::SVec3f, outgoing::SVec3f)::SVec3f
    incoming = reflect(outgoing, normal)
    return incoming
end

@inline function sampleDeltaTransparent(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if (rand(Float32) < fresnelDielectric(ior, normal, outgoing))
        return reflect(outgoing, normal)
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

#####################################
# Eval Delta
#####################################

@inline function evalDelta(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if material.roughness != 0
        return SVec3f(0, 0, 0)
    end

    if material.type == "reflective"
        return evalDeltaReflective(material.color, normal, outgoing)
    elseif material.type == "transparent"
        return evalDeltaTransparent(material.ior, normal, outgoing, incoming)
    elseif material.type == "refractive"
        return evalDeltaRefractive(material.ior, normal, outgoing, incoming)
        # elseif material.type == "volumetric" # not implemented for now
    else
        error("unknown material type")
    end
end

function evalDeltaReflective(
    color::SVec3f,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    radiance =
    #material.color / pi * 
        fresnelConductor(
            reflectivityToEta(color),
            SVec3f(0, 0, 0),
            normal,
            outgoing,
        )# / pi
    return radiance
end

@inline function evalDeltaTransparent(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return vec3f{1,1,1} * fresnelDielectric(ior, normal, outgoing)
    else
        return color * (1 - fresnelDielectric(ior, normal, outgoing))
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
               SVec3f(1, 1, 1) : SVec3f(0, 0, 0)
    end

    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor = entering ? ior : (1.0f0 / ior)
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return SVec3f(1, 1, 1) * fresnelDielectric(relIor, upNormal, outgoing)
    else
        return SVec3f(1, 1, 1) *
               (1.0f0 / (relIor * relIor)) *
               (1 - fresnelDielectric(relIor, upNormal, outgoing))
    end
end

#end module
end