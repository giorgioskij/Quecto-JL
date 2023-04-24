using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra
using ..MaterialsFunctions
using ..Pdfs
using ..Bsdf
using ..Lights

using StaticArrays: dot, cross

export shadeMaterial, shadePathtrace

function shadeMaterial(
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    maxBounces::Int64,
    lights::Vector{Light},
)::SVec3f
    radiance = zeroSV3f
    newRay = ray
    opbounce::Int8 = 0
    weight = SVec3f(1, 1, 1)

    @inbounds for bounce = 1:maxBounces
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance =
                muladd.(
                    weight,
                    evalEnvironment(scene, newRay.direction),
                    radiance,
                )
            break
        end

        instance::Instance = scene.instances[intersection.instanceIndex]
        instanceFrame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]

        outgoing = -newRay.direction
        position = evalShadingPosition(scene, intersection)
        materialPoint::MaterialPoint =
            evalMaterial(scene, instance, intersection)
        # normal = evalShadingNormal(
        #     scene,
        #     instance,
        #     intersection.elementIndex,
        #     intersection.u,
        #     intersection.v,
        #     outgoing,
        # )
        material::Material = scene.materials[instance.materialIndex]

        # normal = evalNormal(
        #     shape,
        #     intersection,
        #     instanceFrame,
        #     outgoing,
        #     materialPoint.type,
        #     material.normalTex,
        # )
        normal = evalNormal(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
            outgoing,
        )

        # handle opacity
        # if !isapprox(materialPoint.opacity, 1) &&
        if materialPoint.opacity < 1 && rand(Float32) >= materialPoint.opacity
            if opbounce > 128
                break
            end
            opbounce += 1
            newRay = Ray(
                muladd.(newRay.direction, 0.01f0, position),
                newRay.direction,
            )
            # newRay = Ray(
            #     position + newRay.direction * 0.01f0,
            #     newRay.direction,
            # )
            bounce -= 1
            continue
        end

        # accumulate emission
        # radiance +=
        #     weight * (
        #         dot(normal, outgoing) >= 0 ? materialPoint.emission :
        #         zeroSV3f
        #     )
        emission =
            dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f
        radiance = muladd.(weight, emission, radiance)

        incoming = zeroSV3f

        # next direction
        if materialPoint.roughness != 0
            # BSDF materials
            incoming = sampleBSDF(materialPoint, normal, outgoing)
            if incoming == zeroSV3f
                break
            end
            weight =
                weight .* evalBSDF(materialPoint, normal, outgoing, incoming) /
                pdfBSDF(materialPoint, normal, outgoing, incoming)
        else
            # Delta materials
            incoming = sampleDelta(materialPoint, normal, outgoing)
            if incoming == zeroSV3f
                break
            end
            weight =
                weight .* evalDelta(materialPoint, normal, outgoing, incoming) /
                pdfDelta(materialPoint, normal, outgoing, incoming)
        end

        # check weight
        if (weight == zeroSV3f || !all(isfinite.(weight)))
            break
        end

        # russian roulette
        if (bounce > 4)
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

function shadePathtrace(
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    maxBounces::Int64,
    lights::Vector{Light},
)::SVec3f
    radiance = zeroSV3f
    newRay = ray
    opbounce::Int8 = 0
    weight = SVec3f(1, 1, 1)

    @inbounds for bounce = 1:maxBounces
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance =
                muladd.(
                    weight,
                    evalEnvironment(scene, newRay.direction),
                    radiance,
                )
            break
        end

        instance::Instance = scene.instances[intersection.instanceIndex]
        instanceFrame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]

        outgoing = -newRay.direction
        position = evalShadingPosition(scene, intersection)
        materialPoint::MaterialPoint =
            evalMaterial(scene, instance, intersection)
        # normal = evalShadingNormal(
        #     scene,
        #     instance,
        #     intersection.elementIndex,
        #     intersection.u,
        #     intersection.v,
        #     outgoing,
        # )
        material::Material = scene.materials[instance.materialIndex]

        # normal = evalNormal(
        #     shape,
        #     intersection,
        #     instanceFrame,
        #     outgoing,
        #     materialPoint.type,
        #     material.normalTex,
        # )
        normal = evalNormal(
            scene,
            instance,
            intersection.elementIndex,
            intersection.u,
            intersection.v,
            outgoing,
        )

        # handle opacity
        # if !isapprox(materialPoint.opacity, 1) &&
        if materialPoint.opacity < 1 && rand(Float32) >= materialPoint.opacity
            if opbounce > 128
                break
            end
            opbounce += 1
            newRay = Ray(
                muladd.(newRay.direction, 0.01f0, position),
                newRay.direction,
            )
            # newRay = Ray(
            #     position + newRay.direction * 0.01f0,
            #     newRay.direction,
            # )
            bounce -= 1
            continue
        end

        # accumulate emission
        # radiance +=
        #     weight * (
        #         dot(normal, outgoing) >= 0 ? materialPoint.emission :
        #         zeroSV3f
        #     )
        emission =
            dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f
        radiance = muladd.(weight, emission, radiance)

        incoming = zeroSV3f

        # next direction
        if materialPoint.roughness != 0
            # BSDF materials
            if rand(Float32) < 0.5f0
                incoming = sampleBSDF(materialPoint, normal, outgoing)
            else
                incoming = sampleLights(scene, lights, position)
            end

            if incoming == zeroSV3f
                break
            end
            weight =
                weight .* evalBSDF(materialPoint, normal, outgoing, incoming) /
                (
                    0.5f0 * pdfBSDF(materialPoint, normal, outgoing, incoming) +
                    0.5f0 * pdfLights(scene, bvh, lights, position, incoming)
                )
        else
            # Delta materials
            incoming = sampleDelta(materialPoint, normal, outgoing)
            # TODO: understand if this if is needed
            if incoming == zeroSV3f
                break
            end
            weight =
                weight .* evalDelta(materialPoint, normal, outgoing, incoming) /
                pdfDelta(materialPoint, normal, outgoing, incoming)
        end

        # check weight
        if (weight == zeroSV3f || !all(isfinite.(weight)))
            break
        end

        # russian roulette
        if (bounce > 4)
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
