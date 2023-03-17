## TODO: to remove if like the new one

using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot, cross

#export shadeMaterial

function shadeMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    radiance = zeroSV3f
    maxBounce = 128
    newRay = ray
    opbounce::Int32 = 0
    weight = SVec3f(1, 1, 1)

    for b = 1:maxBounce
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance += weight * evalEnvironment(scene, newRay.direction)
            break
        end

        instance::Instance = scene.instances[intersection.instanceIndex]
        instanceFrame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]

        outgoing = -newRay.direction
        position = evalShadingPosition(scene, intersection)
        materialPoint::MaterialPoint =
            evalMaterial(scene, instance, intersection)
        normal = evalNormal(
            shape,
            intersection,
            instanceFrame,
            outgoing,
            materialPoint.type,
        )

        # handle opacity
        if !isapprox(materialPoint.opacity, 1) &&
           (rand(Float32) >= materialPoint.opacity)
            if opbounce > 128
                break
            end
            opbounce += 1
            newRay = Ray(position + newRay.direction * 1.0f-2, newRay.direction)
            bounce -= 1
            continue
        end

        # accumulate emission
        radiance +=
            weight *
            (dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f)

        # next direction
        incoming = sampleBSDF(materialPoint, normal, outgoing)
        if incoming == zeroSV3f
            break
        end
        weight *=
            evalBSDF(materialPoint, normal, outgoing, incoming) /
            pdfBSDF(materialPoint, normal, outgoing, incoming)

        # check weight
        if (weight == zeroSV3f || !all(isfinite.(weight)))
            break
        end

        # russian roulette
        if (b > 4)
            rrProb = min(0.99f0, maximum(weight))
            if rand(Float32) >= rrProb
                break
            end
            weight *= (1 / rrProb)
        end

        newRay = Ray(position, incoming)
    end
    return radiance
end

struct MaterialPoint
    type::String
    emission::SVec3f
    color::SVec3f
    opacity::Float32
    roughness::Float32
    ior::Float32

    MaterialPoint(
        type = "matte",
        emission = zeroSV3f,
        color = zeroSV3f,
        opacity = 1,
        roughness = 0,
        ior = 1,
    ) = new(type, emission, color, opacity, roughness, ior)
end

function evalMaterial(
    scene::Scene,
    instance::Instance,
    intersection::Intersection,
)::MaterialPoint
    minRoughness = 0.3f0 * 0.3f0

    material::Material = scene.materials[instance.materialIndex]

    # eval texture coordinates
    textureX, textureY = evalTexcoord(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )

    # eval material textures
    materialEmissionTex::SVec4f =
        evalTexture(scene, material.emissionTex, textureX, textureY)
    materialColorTex::SVec4f =
        evalTexture(scene, material.colorTex, textureX, textureY)
    # ignore roughness and scattering textures

    # eval material properties
    emission::SVec3f = material.emission * xyz(materialEmissionTex)
    color::SVec3f = material.color * xyz(materialColorTex)
    roughness::Float32 = material.roughness

    # fix roughness
    if material.type == "matte" ||
       #material.type == "gltfpbr" ||
       material.type == "glossy"
        roughness = clamp(roughness, minRoughness, 1.0f0)
    elseif roughness < minRoughness
        roughness = 0.0f0
    end

    opacity::Float32 = material.opacity * materialColorTex[4]

    return MaterialPoint(
        material.type,
        emission,
        color,
        opacity,
        roughness,
        material.ior,
    )
end

function pdfBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if material.type == "matte"
        if dot(normal, incoming) * dot(normal, outgoing) <= 0
            return 0
        end
        #upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
        #return pdfHemisphereCos(upNormal, incoming)
        return pdfHemisphereCos(normal, incoming)
    elseif material.type == "reflective"
        if dot(normal, incoming) * dot(normal, outgoing) <= 0
            return 0
        end
        if material.roughness == 0
            return 1.0f0
        else
            #normal = dot(normal, outgoing) <= 0 ? -normal : normal
            halfway = norm(outgoing + incoming)
            return pdfMicrofacet(material.roughness, normal, halfway) /
                   (4.0f0 * abs(dot(outgoing, halfway)))
        end
    elseif material.type == "glossy"
        if dot(normal, incoming) * dot(normal, outgoing) <= 0
            return 0
        end
        #normal = dot(normal, outgoing) <= 0 ? -normal : normal
        F = fresnelDielectric(material.ior, normal, outgoing)
        halfway = norm(outgoing + incoming)

        # if material.roughness == 0
        #     if rand(Float32) < F
        #         return ifelse(
        #             dot(normal, incoming) * dot(normal, outgoing) <= 0,
        #             0.0f0,
        #             1.0f0,
        #         )
        #     else
        #         return F * (4.0f0 * abs(dot(outgoing, halfway))) +
        #                (1 - F) * pdfHemisphereCos(normal, incoming)
        #     end
        # else
        return F * pdfMicrofacet(material.roughness, normal, halfway) /
               (4.0f0 * abs(dot(outgoing, halfway))) +
               (1 - F) * pdfHemisphereCos(normal, incoming)
        #end
    else
        error("Unknown material type")
    end
end

@inline function pdfHemisphereCos(normal::SVec3f, direction::SVec3f)::Float32
    cosw = dot(normal, direction)
    return cosw <= 0 ? 0.0f0 : cosw / pi
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

function sampleBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    # if (material.roughness == 0)
    #     return zeroSV3f
    # end

    if material.type == "matte"
        #upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
        #incoming = sampleHemisphereCos(upNormal)
        incoming = sampleHemisphereCos(normal)

        return incoming
    elseif material.type == "reflective"
        #normal = dot(normal, outgoing) <= 0 ? -normal : normal
        if material.roughness == 0
            incoming = reflect(outgoing, normal)
        else
            # exponent = 2.0f0 / material.roughness^2
            # halfway = sampleHemisphereCosPower(exponent, normal)
            # incoming = reflect(outgoing, halfway)
            halfway = sampleMicrofacet(material.roughness, normal)
            incoming = reflect(outgoing, halfway)
            if (!sameHemisphere(normal, outgoing, incoming))
                return zeroSV3f
            end
        end
        return incoming
    elseif material.type == "glossy"
        #upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
        #F1 = fresnelDielectric(material.ior, upNormal, outgoing)
        F1 = fresnelDielectric(material.ior, normal, outgoing)
        # if material.roughness == 0
        #     if (rand(Float32) < F1)
        #         incoming = reflect(outgoing, normal)
        #         if (!sameHemisphere(normal, outgoing, incoming))
        #             return zeroSV3f
        #         end
        #     else
        #         incoming = sampleHemisphereCos(normal)
        #     end

        #     return incoming
        # else
        if (rand(Float32) < F1)
            #halfway = sampleMicrofacet(material.roughness, upNormal)
            halfway = sampleMicrofacet(material.roughness, normal)
            incoming = reflect(outgoing, halfway)
            #if (!sameHemisphere(upNormal, outgoing, incoming))
            if (!sameHemisphere(normal, outgoing, incoming))
                return zeroSV3f
            end
            return incoming
        else
            # incoming = sampleHemisphereCosPower(
            #     2.0f0 / (material.roughness * material.roughness),
            #     normal,
            # )
            #incoming = sampleHemisphereCos(upNormal)
            incoming = sampleHemisphereCos(normal)
            return incoming
        end
        #end
    else
        error("Unknown material type")
    end
end

function evalBSDF(
    material::MaterialPoint,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return zeroSV3f
    end

    if material.type == "matte"
        #good floor and overall illumination, bad lighting
        radiance = material.color / pi * abs(dot(normal, incoming))

        # good lighting spots, bad floor and overall illumination
        # radiance = linInterp(
        #     radiance, material.color * abs(dot(normal, incoming)),
        #     weight )
        return radiance

    elseif material.type == "reflective"
        #normal = dot(normal, outgoing) <= 0 ? -normal : normal
        if material.roughness == 0
            radiance =
            #material.color / pi * 
                fresnelConductor(
                    reflectivityToEta(material.color),
                    zeroSV3f,
                    normal,
                    outgoing,
                )# / pi

        else
            if dot(normal, incoming) * dot(normal, outgoing) <= 0
                return zeroSV3f
            end
            halfway = norm(incoming + outgoing)
            F = fresnelConductor(
                reflectivityToEta(material.color),
                zeroSV3f,
                halfway,
                incoming,
            )
            D = microfacetDistribution(material.roughness, normal, halfway)
            # D = preciseMicrofacetDistribution(material.roughness, normal, halfway)
            G = microfacetShadowing(
                material.roughness,
                normal,
                halfway,
                outgoing,
                incoming,
            )
            radiance =
            #material.color / pi * 
                F * D * G /
                (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) * #* pi
                abs(dot(normal, incoming))
        end
        return radiance

    elseif material.type == "glossy"
        if dot(normal, incoming) * dot(normal, outgoing) <= 0
            return zeroSV3f
        end
        #normal = dot(normal, outgoing) <= 0 ? -normal : normal
        F1 = fresnelDielectric(material.ior, normal, outgoing)
        halfway = norm(incoming + outgoing)

        # if material.roughness == 0
        #     F = fresnelDielectric(material.ior, halfway, incoming) # 2 * 
        #     D = microfacetDistribution(material.roughness, normal, halfway)
        #     # D = preciseMicrofacetDistribution(material.roughness, normal, halfway)
        #     G = microfacetShadowing(
        #         material.roughness,
        #         normal,
        #         halfway,
        #         outgoing,
        #         incoming,
        #     )

        #     radiance =
        #         material.color * (1 - F1) / (pi) * #2 *  
        #         abs(dot(normal, incoming)) +
        #         SVec3f(1, 1, 1) * F * D * G /
        #         (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
        #         abs(dot(normal, incoming)) #^2.5

        #     return radiance
        # else
        F = fresnelDielectric(material.ior, halfway, incoming)
        D = microfacetDistribution(material.roughness, normal, halfway)
        # D = preciseMicrofacetDistribution(material.roughness, normal, halfway)
        G = microfacetShadowing(
            material.roughness,
            normal,
            halfway,
            outgoing,
            incoming,
        )
        radiance =
            material.color * (1 - F1) / (pi) * #2 *
            abs(dot(normal, incoming)) +
            SVec3f(1, 1, 1) * F * D * G /
            (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
            abs(dot(normal, incoming)) #^1.3
        return radiance
        #end
    else
        error("Unknown material type")
    end
end

@inline function sameHemisphere(
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)
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

@inline function reflect(w::SVec3f, n::SVec3f)::SVec3f
    return -w + 2.0f0 * dot(n, w) * n
end

@inline function fresnelDielectric(
    eta::Float32,
    normal::SVec3f,
    outgoing::SVec3f,
)::Float32
    cosw = abs(dot(normal, outgoing))

    sin2 = 1 - cosw * cosw
    eta2 = eta * eta

    cos2t = 1 - sin2 / eta2
    if cos2t < 0
        return 1
    end

    t0 = sqrt(cos2t)
    t1 = eta * t0
    t2 = eta * cosw

    rs = (cosw - t1) / (cosw + t1)
    rp = (t0 - t2) / (t0 + t2)

    return (rs * rs + rp * rp) / 2
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
    a2plusb2::SVec3f = sqrt.(t0 * t0 + 4.0f0 * eta2 * etak2)
    t1::SVec3f = a2plusb2 .+ cos2
    a::SVec3f = sqrt.((a2plusb2 + t0) / 2.0f0)
    t2::SVec3f = 2.0f0 * a * cosw
    rs::SVec3f = (t1 - t2) / (t1 + t2)

    t3::SVec3f = cos2 * a2plusb2 .+ sin2 * sin2
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
            (cosine2 * roughness2 + 1 - cosine2) *
            (cosine2 * roughness2 + 1 - cosine2)
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
    k = roughness / (dot(nxh, nxh) + a * a)
    return k * k / pi
end

@inline function microfacetShadowing(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
    ggx::Bool = true,
)
    return microfacetShadowing1(roughness, normal, halfway, outgoing, ggx) *
           microfacetShadowing1(roughness, normal, halfway, incoming, ggx)
end

@inline function microfacetShadowing1(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    direction::SVec3f,
    ggx::Bool,
)
    cosine = dot(normal, direction)
    cosineh = dot(halfway, direction)
    if cosine * cosineh <= 0
        return 0.0f0
    end
    roughness2 = roughness * roughness
    cosine2 = cosine * cosine
    if ggx
        return 2 * abs(cosine) /
               (abs(cosine) + sqrt(cosine2 - roughness2 * cosine2 + roughness2))
    else
        ci = abs(cosine) / (roughness * sqrt(1 - cosine2))
        return ci < 1.6f0 ?
               (3.535f0 * ci + 2.181f0 * ci * ci) /
               (1.0f0 + 2.276f0 * ci + 2.577f0 * ci * ci) : 1.0f0
    end
end

@inline function sampleHemisphere(normal::SVec3f)::SVec3f
    ruvx = rand(Float32)
    z = rand(Float32)  # ruvy
    r = sqrt(clamp(1 - z * z, 0.0f0, 1.0f0))
    phi = 2 * pi * ruvx
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

@inline function sampleHemisphereCos(normal::SVec3f)::SVec3f
    ruvx = rand(Float32)
    z = sqrt(rand(Float32))  # ruvy
    r = sqrt(1 - z * z)
    phi = 2 * pi * ruvx
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
    r = sqrt(1.0f0 - z * z)
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
        theta = atan(roughness * sqrt(ruvy / (1 - ruvy)))
    else
        theta = atan(sqrt(-(roughness^2) * log(1 - ruvy)))
    end
    localHalfVector =
        SVec3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
    return transformDirection(basisFromz(normal), localHalfVector)
end

@inline function basisFromz(v::SVec3f)
    z = norm(v)
    sign = copysign(1.0f0, z.z)
    a = -1.0f0 / (sign + z.z)
    b = z.x * z.y * a
    x = SVec3f(1.0f0 + sign * z.x * z.x * a, sign * b, -sign * z.x)
    y = SVec3f(b, sign + z.y * z.y * a, -z.y)
    # return Frame(x, y, z)
    return Mat3f(x, y, z)
end
