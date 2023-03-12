module Pdfs
using ..World
using ..Types
using ..Algebra
using ..MaterialsFunctions

using StaticArrays: dot, cross

export pdfBSDF, pdfDelta

##################################################
# Delta PDFs
##################################################

@inline function pdfDelta(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if material.roughness != 0
        return 0
    end

    if material.type == "reflective"
        return pdfDeltaReflective(normal, outgoing, incoming)
    elseif material.type == "transparent"
        return pdfDeltaTransparent(material.ior, normal, outgoing, incoming)
    elseif material.type == "refractive"
        return pdfDeltaRefractive(material.ior, normal, outgoing, incoming)
        # elseif material.type == "volumetric" # not implemented for now
    else
        error("unknown material type")
    end
end

@inline function pdfDeltaReflective(
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return 0
    else
        return 1.0f0
    end
end

@inline function pdfDeltaTransparent(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return fresnelDielectric(ior, upNormal, outgoing)
    else
        return 1.0f0 - fresnelDielectric(ior, upNormal, outgoing)
    end
end

@inline function pdfDeltaRefractive(
    ior::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if (abs(ior - 1.0f0) < 1.0f-3)
        return dot(normal, incoming) * dot(normal, outgoing) < 0 ? 1.0f0 : 0.0f0
    end
    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor = entering ? ior : (1.0f0 / ior)
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        return fresnelDielectric(relIor, upNormal, outgoing)
    else
        return (1.0f0 - fresnelDielectric(relIor, upNormal, outgoing))
    end
end

################################################
# BSDF PDFs
################################################

@inline function pdfBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if (material.roughness == 0)
        return 0
    end

    if material.type == "matte"
        return pdfBSDFMatte(normal, outgoing, incoming)
    elseif material.type == "glossy"
        pdfBSDFGlossy(
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "reflective"
        pdfBSDFReflective(material.roughness, normal, outgoing, incoming)
    elseif material.type == "transparent"
        pdfBSDFTransparent(
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "refractive"
        pdfBSDFRefractive(
            material.ior,
            material.roughness,
            normal,
            outgoing,
            incoming,
        )
    elseif material.type == "subsurface"
        pdfBSDFRefractive(
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

@inline function pdfBSDFMatte(
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return 0
    end

    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    return pdfHemisphereCos(upNormal, incoming)
end

@inline function pdfBSDFGlossy(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return 0
    end
    F = fresnelDielectric(ior, normal, outgoing)
    halfway = norm(outgoing + incoming)
    return F * pdfMicrofacet(roughness, normal, halfway) /
           (4.0f0 * abs(dot(outgoing, halfway))) +
           (1 - F) * pdfHemisphereCos(normal, incoming)
end

@inline function pdfBSDFReflective(
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return 0
    end

    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    halfway = norm(outgoing + incoming)
    return pdfMicrofacet(roughness, upNormal, halfway) /
           (4.0f0 * abs(dot(outgoing, halfway)))
end

@inline function pdfBSDFTransparent(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        halfway = norm(incoming + outgoing)
        return fresnelDielectric(ior, halfway, outgoing) *
               pdfMicrofacet(roughness, upNormal, halfway) /
               (4.0f0 * abs(dot(outgoing, halfway)))
    else
        reflected = reflect(-incoming, upNormal)
        halfway = norm(reflected + outgoing)
        d =
            (1 - fresnelDielectric(ior, halfway, outgoing)) *
            pdfMicrofacet(roughness, upNormal, halfway)
        return d / (4.0f0 * abs(dot(outgoing, halfway)))
    end
end

@inline function pdfBSDFRefractive(
    ior::Float32,
    roughness::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    entering = dot(normal, outgoing) >= 0
    upNormal = entering ? normal : -normal
    relIor = entering ? ior : (1.0f0 / ior)
    if (dot(normal, incoming) * dot(normal, outgoing) >= 0)
        halfway = norm(incoming + outgoing)
        return fresnelDielectric(relIor, halfway, outgoing) *
               pdfMicrofacet(roughness, upNormal, halfway) /
               (4.0f0 * abs(dot(outgoing, halfway)))
    else
        halfway =
            -norm(relIor * incoming + outgoing) * (entering ? 1.0f0 : -1.0f0)
        return (1.0f0 - fresnelDielectric(relIor, halfway, outgoing)) *
               pdfMicrofacet(roughness, upNormal, halfway) *
               abs(dot(halfway, incoming)) /
               (relIor * dot(halfway, incoming) + dot(halfway, outgoing))^2.0f0
    end
end

@inline function pdfHemisphereCos(normal::SVec3f, direction::SVec3f)::Float32
    cosw = dot(normal, direction)
    return (cosw <= 0) ? 0.0f0 : (cosw / pi)
end

@inline function pdfMicrofacet(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    ggx::Bool = true,
)::Float32
    cosine = dot(normal, halfway)
    if cosine < 0
        return 0
    end
    return microfacetDistribution(roughness, normal, halfway, ggx) * cosine
    # preciseMicrofacetDistribution(roughness, normal, halfway, ggx) * cosine
end

#end module
end