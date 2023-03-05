
using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot, cross

export shadeRaytrace

function shadeRaytrace(
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