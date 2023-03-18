using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra
using ..MaterialsFunctions
using ..Pdfs
using ..Bsdf

using StaticArrays: dot, cross

export shadeMaterial

function shadeMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    radiance = zeroSV3f
    maxBounce = 128
    newRay = ray
    opbounce::Int8 = 0
    weight = SVec3f(1, 1, 1)

    @inbounds for bounce = 1:maxBounce
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance =
                muladd.(
                    weight,
                    evalEnvironment(scene, newRay.direction),
                    radiance,
                )
            #radiance += weight * evalEnvironment(scene, newRay.direction)
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
            newRay = Ray(
                muladd.(newRay.direction, 0.01f0, position),
                newRay.direction,
            )
            #newRay = Ray(position + newRay.direction * 0.01f0, newRay.direction)
            bounce -= 1
            continue
        end

        # accumulate emission
        radiance +=
            weight *
            (dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f)

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
