module Shading

using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot, cross

export shaderEyelight,
    shaderMaterial,
    shaderNormal,
    shaderIndirectNaive,
    shaderEyelightBsdf,
    shaderPath

# function shaderColor(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
#     intersection::Intersection = intersectScene(ray, scene)

#     if !intersection.hit
#         radiance = evalEnvironment(scene, ray.direction)
#         return radiance
#     end

#     instance::Instance = scene.instances[intersection.instanceIndex]
#     material::Material = scene.materials[instance.materialIndex]

#     radiance = material.color

#     return radiance
# end

function shaderNormal(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]

    outgoing = -ray.direction
    normal = evalNormal(shape, intersection, frame)
    if dot(normal, outgoing) < 0
        normal = -normal
    end

    radiance::SVec3f = normal * 0.5 .+ 0.5

    return srgbToRgb(radiance)
end

# this is meant to be a copy of trace_naive but with only matte
function shaderPath(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    radiance::SVec3f = SVec3f(0, 0, 0)
    weight::SVec3f = SVec3f(1, 1, 1)
    hit::Bool = false
    opbounce::Int = 0

    nbounces = 8
    # trace path
    for bounce = 1:nbounces

        # intersect next point
        intersection::Intersection = intersectScene(ray, scene, bvh)
        if !intersection.hit
            # FIX: fuck , in yocto evalEnvironment yields a color with no interval
            radiance += weight * evalEnvironment(scene, ray.direction)
            # println("Bounce $bounce, radiance: $radiance")
            break
        end

        # prepare shading point
        outgoing = -ray.direction
        instance::Instance = scene.instances[intersection.instanceIndex]
        shape::Shape = scene.shapes[instance.shapeIndex]
        position = evalPosition(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
        )
        normal = evalNormal(shape, intersection, instance.frame)

        # no evalMaterial because we do so explicitly
        textureX, textureY = evalTexcoord(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
        )
        material::Material = scene.materials[instance.materialIndex]
        materialEmissionTex::SVec4f =
            evalTexture(scene, material.emissionTex, textureX, textureY)
        materialEmission::SVec3f = material.emission * xyz(materialEmissionTex)

        materialColorTex::SVec4f =
            evalTexture(scene, material.colorTex, textureX, textureY)
        materialColor::SVec3f = material.color * xyz(materialColorTex)

        # roughness texture not implemented for now
        minRoughness = 0.03f0 * 0.03f0
        materialRoughness = material.roughness
        materialRoughness = materialRoughness * materialRoughness
        # fix roughness - only matte for now
        materialRoughness = clamp(materialRoughness, minRoughness, 1.0f0)

        # handle opacity (no need for now, keep it simple)

        # set hit variables
        if bounce == 1
            hit = true
            # don't care about hit albedo, hit normal
        end

        # accumulate emission
        radiance += weight * evalEmission(materialEmission, normal, outgoing)

        # next direction
        incoming::SVec3f = SVec3f(0, 0, 0)
        if materialRoughness != 0
            incoming = sampleBsdfCos(
                material.type,
                materialRoughness,
                materialColor,
                normal,
                outgoing,
            )
            if incoming == SVec3f(0, 0, 0)
                break
            end
            weight *=
                evalBsdfCos(
                    material.type,
                    materialRoughness,
                    materialColor,
                    normal,
                    outgoing,
                    incoming,
                ) / sampleBsdfCosPdf(
                    material.type,
                    materialRoughness,
                    materialColor,
                    normal,
                    outgoing,
                    incoming,
                )
        else
            error("Roughness is not 0")
        end

        # check weight
        if weight == SVec3f(0, 0, 0) || !all(isfinite.(weight))
            break
        end

        # russian roulette
        if bounce > 4
            rrProb = min(0.99f0, maximum(weight))
            if rand(Float32) >= rrProb
                break
            end
            weight *= 1 / rrProb
        end

        # setup next iteration
        ray = Ray(position, incoming)
    end

    return radiance
end

function sampleBsdfCos(
    materialType::String,
    materialRoughness::Float32,
    materialColor::SVec3f,
    normal,
    outgoing,
)::SVec3f
    if materialRoughness == 0
        error("Wait shouldnt be 0")
        return SVec3f(0, 0, 0)
    end
    # only matte 
    if materialType == "matte"
        return sampleMatte(normal, outgoing)
    else
        error("only matte for now")
    end
end
function sampleMatte(normal::SVec3f, outgoing::SVec3f)
    upNormal = dot(normal, outgoing) <= 0 ? -normal : normal
    return sampleHemisphereCos(upNormal)
end
function sampleBsdfCosPdf(
    materialType::String,
    materialRoughness::Float32,
    materialColor::SVec3f,
    normal::SVec3f,
    outgoing::SVec3f,
    incoming::SVec3f,
)::Float32
    if materialRoughness == 0
        error("Wait shouldnt be 0")
        return SVec3f(0, 0, 0)
    end
    # only matte 
    if materialType == "matte"
        return sampleMattePdf(normal, outgoing, incoming)
    else
        error("only matte for now")
    end
end
function sampleMattePdf(normal, ougoing, incoming)::Float32
    if dot(normal, incoming) * dot(normal, ougoing) <= 0
        return 0
    end
    upNormal = dot(normal, ougoing) <= 0 ? -normal : normal
    return sampleHemisphereCosPdf(upNormal, incoming)
end
function sampleHemisphereCosPdf(normal::SVec3f, direction::SVec3f)::Float32
    cosw = dot(normal, direction)
    return cosw <= 0 ? 0 : cosw / pi
end
function evalBsdfCos(
    materialType::String,
    materialRoughness::Float32,
    materialColor::SVec3f,
    normal,
    outgoing,
    incoming,
)::SVec3f
    if materialRoughness == 0
        error("wait what")
        return SVec3f(0, 0, 0)
    end

    # only matte materials
    if materialType == "matte"
        return evalMatte(materialColor, normal, outgoing, incoming)
    else
        error("only matte for now")
    end
end

@inline function evalMatte(materialColor, normal, outgoing, incoming)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return SVec3f(0, 0, 0)
    end
    return materialColor / pi * abs(dot(normal, incoming))
end

function evalEmission(materialEmission, normal, outgoing)::SVec3f
    return dot(normal, outgoing) >= 0 ? materialEmission : SVec3f(0, 0, 0)
end

function shaderEyelightBsdf(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f

    # initialize
    radiance::SVec3f = SVec3f(0, 0, 0)
    weight::SVec3f = SVec3f(1, 1, 1)
    maxBounces = 4
    opacityBounce = 0

    for bounce = 1:maxBounces

        # intersect 
        intersection::Intersection = intersectScene(ray, scene, bvh, false)
        if !intersection.hit
            radiance += weight * evalEnvironment(scene, ray.direction)
            break
        end

        # extract intersection data
        instance::Instance = scene.instances[intersection.instanceIndex]
        shape::Shape = scene.shapes[instance.shapeIndex]
        frame::Frame = instance.frame

        # eval position
        outgoing::SVec3f = -ray.direction
        position::SVec3f = evalPosition(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
        )

        # compute normal with correction
        normal::SVec3f = evalNormal(shape, intersection, frame)
        if dot(normal, outgoing) < 0
            normal = -normal
        end

        # compute texture coordinates
        textureX, textureY = evalTexcoord(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
        )

        # eval color of the shape - by default {1,1,1,1}
        shapeColor::SVec4f = SVec4f(1, 1, 1, 1)

        # evaluate color and emission of the  material
        material::Material = scene.materials[instance.materialIndex]

        # emission
        materialEmissionTex::SVec4f =
            evalTexture(scene, material.emissionTex, textureX, textureY)
        materialEmission::SVec3f =
            material.emission * xyz(materialEmissionTex) * xyz(shapeColor)

        # color
        materialColorTex::SVec4f =
            evalTexture(scene, material.colorTex, textureX, textureY)
        materialColor::SVec3f =
            material.color * xyz(materialColorTex) * xyz(shapeColor)

        # evaluate opacity
        materialOpacity = material.opacity * materialColorTex[4] * shapeColor[4]

        # fix minimum roughness
        minRoughness = 0.03f0 * 0.03f0
        materialRoughness = material.roughness
        # if material.type == "matte"
        #     error("its matte")
        materialRoughness = clamp(material.roughness, minRoughness, 1.0f0)
        # end

        # missing: metallic,  ior, scattering...

        # handle opacity
        if !isapprox(materialOpacity, 1) && rand(Float32) >= materialOpacity
            if opacityBounce > 128
                break
            end
            opacityBounce += 1
            ray = Ray(position + ray.direction * 1.0f-2, ray.direction)
            bounce -= 1
            continue
        end

        # set hit variables
        if bounce == 1
            hit = true
            # missing: hit albedo, hit normal
        end

        # accumulate emission
        incoming::SVec3f = outgoing

        # missing: bsdf
        emission =
            dot(normal, outgoing) >= 0 ? materialEmission : SVec3f(0, 0, 0)
        radiance += weight * emission

        # brdf + light
        radiance +=
            (weight * pi) * evalBsdfCos(
                material.type,
                materialRoughness,
                materialColor,
                normal,
                outgoing,
                incoming,
            )

        # if material is not delta, break
        break
    end
    return radiance
end

function shaderEyelight(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # EVALUATE GEOMETRY 
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]
    normal = evalNormal(shape, intersection, frame)

    textureX, textureY = evalTexcoord(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )
    outgoing = -ray.direction

    # EVALUATE MATERIAL
    material::Material = scene.materials[instance.materialIndex]
    emission::SVec3f =
        material.emission *
        xyz(evalTexture(scene, material.emissionTex, textureX, textureY))
    color::SVec3f =
        material.color *
        xyz(evalTexture(scene, material.colorTex, textureX, textureY))

    # materialColor = evalMaterialColor(scene, intersection)

    # radiance::SVec3f = emission + abs(dot(normal, outgoing)) * color
    radiance::SVec3f = abs(dot(normal, outgoing)) * color

    return SVec3f(radiance.x, radiance.y, radiance.z)
end

function shaderIndirectNaive(
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    bounce::Int = 1,
)::SVec3f
    maxBounce = 32

    # INTERSECT SCENE
    intersection::Intersection = intersectScene(ray, scene, bvh)

    # EVAL ENVIRONMENT
    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # EVALUATE GEOMETRY
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]

    # position = evalShadingPosition(scene, intersection, outgoing)
    position = evalPosition(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )
    # normal = evalNormal(shape, intersection, frame)
    # FIX: important! here it's instance.frame, not intersection.frame !
    normal = evalNormal(shape, intersection, instance.frame)
    textureX, textureY = evalTexcoord(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )
    outgoing = -ray.direction

    # NORMAL CORRECTIONS
    if dot(normal, outgoing) < 0
        normal = -normal
    end

    # EVALUATE MATERIAL
    material::Material = scene.materials[instance.materialIndex]
    materialEmissionTex::SVec4f =
        evalTexture(scene, material.emissionTex, textureX, textureY)
    materialColorTex::SVec4f =
        evalTexture(scene, material.colorTex, textureX, textureY)
    emission::SVec3f = material.emission * xyz(materialEmissionTex)
    color::SVec3f = material.color * xyz(materialColorTex)
    roughness = material.roughness * material.roughness
    opacity::Float32 = material.opacity * materialColorTex[4]

    # HANDLE OPACITY
    if rand(Float32) >= opacity
        return shaderIndirectNaive(
            scene,
            Ray(position, ray.direction),
            bvh,
            bounce + 1,
        )
    end

    # ACCUMULATE EMISSION
    radiance = emission

    # EXIT IF RAY IS DONE
    if bounce >= maxBounce
        return radiance
    end

    # COMPUTE ILLUMINATION 

    if material.type == "matte"
        upNormal = ifelse(dot(normal, outgoing) <= 0, -normal, normal)
        incoming = sampleHemisphereCos(normal)
        if incoming == SVec3f(0, 0, 0)
            return radiance
        end

        #good floor and overall illumination, bad lighting
        radiance += ifelse(
            dot(upNormal, incoming) * dot(upNormal, outgoing) <= 0,
            SVec3f(0, 0, 0),
            color *
            abs(dot(upNormal, incoming)) *
            shaderIndirectNaive(scene, Ray(position, incoming), bvh, bounce + 1),
        )
        # radiance +=
        #     color *
        #     shaderIndirectNaive(scene, Ray(position, incoming), bvh, bounce + 1)

        # good lighting spots, bad floor and overall illumination
        # radiance = linInterp(
        #     radiance,
        #     ifelse(
        #         dot(normal, incoming) * dot(normal, outgoing) <= 0,
        #         SVec3f(0, 0, 0),
        #         color * abs(dot(normal, incoming)),
        #     ),
        #     weight,
        # )

    elseif material.type == "reflective"
        if material.roughness == 0
            upNormal = ifelse(dot(normal, outgoing) <= 0, -normal, normal)
            incoming = reflect(outgoing, upNormal)
            if incoming == SVec3f(0, 0, 0)
                return radiance
            end

            radiance += ifelse(
                dot(upNormal, incoming) * dot(upNormal, outgoing) <= 0,
                SVec3f(0, 0, 0),
                # color *
                fresnelConductor(
                    reflectivityToEta(color),
                    SVec3f(0, 0, 0),
                    upNormal,
                    outgoing,
                ) * shaderIndirectNaive(
                    scene,
                    Ray(position, incoming),
                    bvh,
                    bounce + 1,
                ),
            )

        else
            upNormal = ifelse(dot(normal, outgoing) <= 0, -normal, normal)
            halfway = sampleMicrofacet(material.roughness, upNormal)
            incoming = reflect(outgoing, halfway)
            if (!sameHemisphere(upNormal, outgoing, incoming))
                return radiance
            end
            if incoming == SVec3f(0, 0, 0)
                return radiance
            end
            #halfway = norm(incoming + outgoing)
            F = fresnelConductor(
                reflectivityToEta(color),
                SVec3f(0, 0, 0),
                halfway,
                incoming,
            )
            D = microfacetDistribution(material.roughness, upNormal, halfway)
            G = microfacetShadowing(
                material.roughness,
                upNormal,
                halfway,
                outgoing,
                incoming,
            )
            radiance += ifelse(
                dot(upNormal, incoming) * dot(upNormal, outgoing) <= 0,
                SVec3f(0, 0, 0),
                #color * 
                F * D * G / (
                    4.0f0 * dot(upNormal, outgoing) * dot(upNormal, incoming)
                ) *
                abs(dot(upNormal, incoming)) *
                shaderIndirectNaive(
                    scene,
                    Ray(position, incoming),
                    bvh,
                    bounce + 1,
                ),
            )
        end
    elseif material.type == "glossy"
        upNormal = ifelse(dot(normal, outgoing) <= 0, -normal, normal)
        F1 = fresnelDielectric(material.ior, upNormal, outgoing)
        if material.roughness == 0
            if (rand(Float32) < F1)
                incoming = reflect(outgoing, upNormal)
                if (!sameHemisphere(upNormal, outgoing, incoming))
                    return radiance
                end
            else
                incoming = sampleHemisphereCos(upNormal)
            end
            if incoming == SVec3f(0, 0, 0)
                return radiance
            end

            halfway = norm(incoming + outgoing)
            F = fresnelDielectric(2 * material.ior, halfway, incoming)
            # how we compute D for sharp plastic?
            #D = microfacetDistribution(material.roughness, upNormal, halfway)
            G = microfacetShadowing(
                material.roughness,
                upNormal,
                halfway,
                outgoing,
                incoming,
            )
            radiance += ifelse(
                dot(upNormal, incoming) * dot(upNormal, outgoing) <= 0,
                SVec3f(0, 0, 0),
                color * (1 - F1) / pi * abs(dot(upNormal, incoming)) +
                SVec3f(1, 1, 1) * F * G / (
                    4.0f0 * dot(upNormal, outgoing) * dot(upNormal, incoming)
                ) *
                abs(dot(upNormal, incoming))^1.75 *
                shaderIndirectNaive(
                    scene,
                    Ray(position, incoming),
                    bvh,
                    bounce + 1,
                ),
            )

        else
            if (rand(Float32) < F1)
                halfway = sampleMicrofacet(material.roughness, upNormal)
                incoming = reflect(outgoing, halfway)
                if (!sameHemisphere(upNormal, outgoing, incoming))
                    return radiance
                end
            else
                incoming = sampleHemisphereCosPower(
                    2.0f0 / (material.roughness * material.roughness),
                    upNormal,
                )
                #incoming = sampleHemisphereCos(upNormal)
            end

            if incoming == SVec3f(0, 0, 0)
                return radiance
            end

            halfway = norm(incoming + outgoing)
            F = fresnelDielectric(material.ior, halfway, incoming)
            D = microfacetDistribution(material.roughness, upNormal, halfway)
            #D = preciseMicrofacetDistribution(material.roughness, upNormal, halfway)
            G = microfacetShadowing(
                material.roughness,
                upNormal,
                halfway,
                outgoing,
                incoming,
            )
            radiance += ifelse(
                dot(upNormal, incoming) * dot(upNormal, outgoing) <= 0,
                SVec3f(0, 0, 0),
                color * (1 - F1) / pi * abs(dot(upNormal, incoming)) +
                SVec3f(1, 1, 1) * F * D * G / (
                    4.0f0 * dot(upNormal, outgoing) * dot(upNormal, incoming)
                ) *
                abs(dot(upNormal, incoming))^1.05 *
                shaderIndirectNaive(
                    scene,
                    Ray(position, incoming),
                    bvh,
                    bounce + 1,
                ),
            )
        end

    else
        error("Unknown material type")
    end
    return radiance
end

function shaderMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    radiance = SVec3f(0, 0, 0)
    maxBounce = 128
    newRay = ray
    opbounce::Int32 = 0
    weight = SVec3f(1, 1, 1)

    for b = 1:maxBounce
        # weight::Float32 = 1.0f0 / b
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            #if b > 1
            radiance += weight * evalEnvironment(scene, newRay.direction)
            #end
            break
        end

        # compute normal of the point hit
        instance::Instance = scene.instances[intersection.instanceIndex]
        frame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]
        material::Material = scene.materials[instance.materialIndex]
        textureX, textureY = evalTexcoord(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
        )

        outgoing = -newRay.direction
        position = evalShadingPosition(scene, intersection, outgoing)
        normal = evalNormal(shape, intersection, frame)
        if dot(normal, outgoing) <= 0
            normal = -normal
        end

        materialEmissionTex::SVec4f =
            evalTexture(scene, material.emissionTex, textureX, textureY)
        materialColorTex::SVec4f =
            evalTexture(scene, material.colorTex, textureX, textureY)
        emission::SVec3f = material.emission * xyz(materialEmissionTex)
        color::SVec3f = material.color * xyz(materialColorTex)

        #incoming = outgoing   # decomment for eyelight shader
        opacity::Float32 = material.opacity * materialColorTex[4]

        # HANDLE OPACITY
        if (material.opacity < 1 && rand(Float32) >= opacity)
            if opbounce > 128
                break
            end
            opbounce += 1
            newRay = Ray(position + newRay.direction * 1e-2f, newRay.direction)
            bounce -= 1
            continue
        end

        #radiance += weight * emission

        # if (rand(Float32) < 0.5f0) # TODO: implement lights
        incoming = sampleBSDF(material, normal, outgoing)
        # else
        #     incoming = sample_lights(scene, lights, position) # TODO: implement lights
        # end

        if incoming == SVec3f(0, 0, 0)
            break
        end

        weight *=
            evalBSDF(material, normal, incoming, outgoing, color) /
            pdfBSDF(material, normal, outgoing, incoming)

        radiance =
            weight * evalBSDF(material, normal, incoming, outgoing, color)

        newRay = Ray(position, incoming)
    end
    return radiance
end

function pdfBSDF(material, normal, outgoing, incoming)
    if material.type == "matte"
        return pdfHemisphereCos(normal, incoming)
    elseif material.type == "reflective"
        if material.roughness == 0
            return ifelse(
                dot(normal, incoming) * dot(normal, outgoing) <= 0,
                0.0f0,
                1.0f0,
            )
        else
            halfway = norm(incoming + outgoing)
            return pdfMicrofacet(material.roughness, normal, halfway) /
                   (4.0f0 * abs(dot(outgoing, halfway)))
        end
    elseif material.type == "glossy"
        F = fresnelDielectric(material.ior, normal, outgoing)
        halfway = norm(incoming + outgoing)

        if material.roughness == 0
            if rand(Float32) < F
                return ifelse(
                    dot(normal, incoming) * dot(normal, outgoing) <= 0,
                    0.0f0,
                    1.0f0,
                )
            else
                return F * (4.0f0 * abs(dot(outgoing, halfway))) +
                       (1 - F) * pdfHemisphereCos(normal, incoming)
            end
        else
            return F * pdfMicrofacet(material.roughness, normal, halfway) /
                   (4.0f0 * abs(dot(outgoing, halfway))) +
                   (1 - F) * pdfHemisphereCos(normal, incoming)
        end
    else
        error("Unknown material type")
    end
end

@inline function pdfHemisphereCos(normal::SVec3f, direction::SVec3f)::Float32
    cosw = dot(normal, direction)
    return ifelse(cosw <= 0, 0.0f0, cosw / pi)
end

@inline function pdfMicrofacet(
    roughness::Float32,
    normal::SVec3f,
    halfway::SVec3f,
    ggx::Bool = true,
)::Float32
    cosine = dot(normal, halfway)
    return ifelse(
        cosine <= 0,
        0.0f0,
        # preciseMicrofacetDistribution(roughness, normal, halfway, ggx) * cosine
        microfacetDistribution(roughness, normal, halfway, ggx) * cosine,
    )
end

function sampleBSDF(
    material::Material,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    # if (material.roughness == 0)
    #     return SVec3f(0, 0, 0)
    # end

    if material.type == "matte"
        incoming = sampleHemisphereCos(normal)

        return incoming
    elseif material.type == "reflective"
        if material.roughness == 0
            incoming = reflect(outgoing, normal)
        else
            # exponent = 2.0f0 / material.roughness^2
            # halfway = sampleHemisphereCosPower(exponent, normal)
            # incoming = reflect(outgoing, halfway)
            halfway = sampleMicrofacet(material.roughness, normal)
            incoming = reflect(outgoing, halfway)
            if (!sameHemisphere(normal, outgoing, incoming))
                return SVec3f(0, 0, 0)
            end
        end
        return incoming
    elseif material.type == "glossy"
        F1 = fresnelDielectric(material.ior, normal, outgoing)
        if material.roughness == 0
            if (rand(Float32) < F1)
                incoming = reflect(outgoing, normal)
                if (!sameHemisphere(normal, outgoing, incoming))
                    return SVec3f(0, 0, 0)
                end
            else
                incoming = sampleHemisphereCos(normal)
            end

            return incoming
        else
            if (rand(Float32) < F1)
                halfway = sampleMicrofacet(material.roughness, normal)
                incoming = reflect(outgoing, halfway)
                if (!sameHemisphere(normal, outgoing, incoming))
                    return SVec3f(0, 0, 0)
                end
            else
                incoming = sampleHemisphereCosPower(
                    2.0f0 / (material.roughness * material.roughness),
                    normal,
                )
                #incoming = sampleHemisphereCos(upNormal)
            end

            return incoming
        end
    else
        error("Unknown material type")
    end
end

function evalBSDF(material, normal, incoming, outgoing, color)
    if dot(normal, incoming) * dot(normal, outgoing) <= 0
        return SVec3f(0, 0, 0)
    end

    if material.type == "matte"
        #good floor and overall illumination, bad lighting
        radiance = color / pi * abs(dot(normal, incoming))

        # good lighting spots, bad floor and overall illumination
        # radiance = linInterp(
        #     radiance, color * abs(dot(normal, incoming)),
        #     weight )
        return radiance

    elseif material.type == "reflective"
        if material.roughness == 0
            radiance =
            #color / pi * 
                fresnelConductor(
                    reflectivityToEta(color),
                    SVec3f(0, 0, 0),
                    normal,
                    outgoing,
                ) / pi

        else
            halfway = norm(incoming + outgoing)
            F = fresnelConductor(
                reflectivityToEta(color),
                SVec3f(0, 0, 0),
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
            #color / pi * 
                F * D * G /
                (4.0f0 * pi * dot(normal, outgoing) * dot(normal, incoming)) *
                abs(dot(normal, incoming))
        end
        return radiance

    elseif material.type == "glossy"
        F1 = fresnelDielectric(material.ior, normal, outgoing)
        halfway = norm(incoming + outgoing)

        if material.roughness == 0
            F = fresnelDielectric(2 * material.ior, halfway, incoming)
            # D = microfacetDistribution(material.roughness, normal, halfway)
            # D = preciseMicrofacetDistribution(material.roughness, normal, halfway)
            # TODO: how to compute D for sharp plastic?
            G = microfacetShadowing(
                material.roughness,
                normal,
                halfway,
                outgoing,
                incoming,
            )

            radiance =
                color * (1 - F1) / (2 * pi) * abs(dot(normal, incoming)) +
                SVec3f(1, 1, 1) * F * G /
                (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
                abs(dot(normal, incoming))^2.5

            return radiance
        else
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
                color * (1 - F1) / (2 * pi) * abs(dot(normal, incoming)) +
                SVec3f(1, 1, 1) * F * D * G /
                (4.0f0 * dot(normal, outgoing) * dot(normal, incoming)) *
                abs(dot(normal, incoming))^1.3
            return radiance
        end
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
#     if specular == SVec3f(0, 0, 0)
#         return SVec3f(0, 0, 0)
#     end
#     cosine = dot(normal, outgoing)
#     specular .+
#     (1.0f0 .- specular) .* clamp(1.0f0 .- abs(cosine), 0.0f0, 1.0f0)^5.0f0
# end

@inline function reflectivityToEta(color::SVec3f)::SVec3f
    reflectivity::SVec3f = clamp.(color, 0.0f0, 1.0f0)
    return (1.0f0 .+ sqrt.(reflectivity)) ./ (1.0f0 .- sqrt.(reflectivity))
end

@inline function reflect(w::SVec3f, n::SVec3f)::SVec3f
    return -w + 2.0f0 * dot(w, n) * n
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
        return SVec3f(0, 0, 0)
    end
    cosw = clamp(cosw, -1.0f0, 1.0f0)
    cos2::Float32 = cosw * cosw
    sin2::Float32 = clamp(1 - cos2, 0.0f0, 1.0f0)
    eta2::SVec3f = eta .* eta
    etak2::SVec3f = etak .* etak

    t0::SVec3f = eta2 .- etak2 .- sin2
    a2plusb2::SVec3f = sqrt.(t0 .* t0 .+ 4.0f0 .* eta2 .* etak2)
    t1::SVec3f = a2plusb2 .+ cos2
    a::SVec3f = sqrt.((a2plusb2 .+ t0) ./ 2.0f0)
    t2::SVec3f = 2.0f0 .* a .* cosw
    rs::SVec3f = (t1 .- t2) ./ (t1 .+ t2)

    t3::SVec3f = cos2 .* a2plusb2 .+ sin2 .* sin2
    t4::SVec3f = t2 .* sin2
    rp::SVec3f = rs .* (t3 .- t4) ./ (t3 .+ t4)

    return (rp .+ rs) ./ 2.0f0
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
        return roughness2 / ((
            pi *
            (cosine2 * roughness2 + 1 - cosine2) *
            (cosine2 * roughness2 + 1 - cosine2)
        ))
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
        return ifelse(
            ci < 1.6f0,
            (3.535f0 * ci + 2.181f0 * ci * ci) /
            (1.0f0 + 2.276f0 * ci + 2.577f0 * ci * ci),
            1.0f0,
        )
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
        theta = atan(sqrt(-roughness^2 * log(1 - ruvy)))
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
    return Frame(x, y, z)
end

end