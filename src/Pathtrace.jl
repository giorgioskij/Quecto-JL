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

export shadeMaterial, shadePath, shadeVolumetric

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
        radiance +=
            weight *
            (dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f)
        # emission =
        #     dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f
        # radiance = muladd.(weight, emission, radiance)

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

function shadePath(
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
            # radiance =
            # muladd.(
            #     weight,
            #     evalEnvironment(scene, newRay.direction),
            #     radiance,
            # )
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
        radiance +=
            weight *
            (dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f)
        # emission =
        #     dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f
        # radiance = muladd.(weight, emission, radiance)

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

@inline function sampleTransmittance(
    density::SVec3f,
    maxDistance::Float32,
)::Float32
    channel = clamp(floor(Int32, rand(Float32) * 3), 0, 2) + 1
    distance = (
        density[channel] == 0 ? Base.max_values(Float32) :
        -log(1 - rand(Float32) / density[channel])
    )
    return min(distance, maxDistance)
end

@inline function sampleTransmittancePdf(
    density::SVec3f,
    distance::Float32,
    maxDistance::Float32,
)::Float32
    if (distance < maxDistance)
        return sum(density * exp.(-density * distance)) / 3
    else
        return sum(exp.(-density * maxDistance)) / 3
    end
end

@inline function evalTransmittance(density::SVec3f, distance::Float32)
    return exp.(-density * distance)
end

function evalScattering(
    material::MaterialPoint,
    outgoing::SVec3f,
    incoming::SVec3f,
)::SVec3f
    if material.density == zeroSV3f
        return zeroSV3f
    end
    return material.scattering *
           material.density *
           evalPhasefunction(material.scanisotropy, outgoing, incoming)
end

function sampleScattering(material::MaterialPoint, outgoing::SVec3f)::SVec3f
    if material.density == zeroSV3f
        return zeroSV3f
    end
    return samplePhasefunction(material.scanisotropy, outgoing)
end

function sampleScatteringPdf(
    material::MaterialPoint,
    outgoing::SVec3f,
    incoming::SVec3f,
)::Float32
    if material.density == zeroSV3f
        return 0
    end
    return samplePhasefunctionPdf(material.scanisotropy, outgoing, incoming)
end

@inline function evalPhasefunction(
    anisotropy::Float32,
    outgoing::SVec3f,
    incoming::SVec3f,
)::Float32
    cosine = -dot(outgoing, incoming)
    denom = 1 + anisotropy * anisotropy - 2 * anisotropy * cosine
    return (1 - anisotropy * anisotropy) / (4 * pi * denom * sqrt(denom))
end

@inline function samplePhasefunction(
    anisotropy::Float32,
    outgoing::SVec3f,
)::SVec3f
    cosTheta = 0.0f0
    rnx = rand(Float32)
    rny = rand(Float32)
    if abs(anisotropy) < 1.0f-3
        cosTheta = 1 - 2 * rny
    else
        square =
            (1 - anisotropy * anisotropy) /
            (1 + anisotropy - 2 * anisotropy * rny)
        cosTheta =
            (1 + anisotropy * anisotropy - square * square) / (2 * anisotropy)
    end

    sinTheta = sqrt(max(0.0f0, 1 - cosTheta * cosTheta))
    phi = 2 * pi * rnx
    localIncoming = SVec3f(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta)
    return matMulVec(basisFromz(-outgoing), localIncoming)
end

@inline function samplePhasefunctionPdf(
    anisotropy::Float32,
    outgoing::SVec3f,
    incoming::SVec3f,
)::Float32
    return evalPhasefunction(anisotropy, outgoing, incoming)
end

function shadeVolumetric(
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
    volumeStack = Vector{MaterialPoint}(undef, 0)

    @inbounds for bounce = 1:maxBounces
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance += weight * evalEnvironment(scene, newRay.direction)
            break
        end

        # handle transmission if inside a volume
        inVolume::Bool = false
        if (!isempty(volumeStack))
            vsdf = volumeStack[end]
            distance = sampleTransmittance(vsdf.density, intersection.distance)
            weight *=
                evalTransmittance(vsdf.density, distance) /
                sampleTransmittancePdf(
                    vsdf.density,
                    distance,
                    intersection.distance,
                )
            inVolume = distance < intersection.distance
            intersection = Intersection(
                intersection.hit,
                intersection.instanceIndex,
                intersection.elementIndex,
                intersection.u,
                intersection.v,
                distance,
            )
        end

        if !inVolume
            instance::Instance = scene.instances[intersection.instanceIndex]
            outgoing = -newRay.direction
            position = evalShadingPosition(scene, intersection)
            materialPoint::MaterialPoint =
                evalMaterial(scene, instance, intersection)

            normal = evalNormal(
                scene,
                instance,
                intersection.elementIndex,
                intersection.u,
                intersection.v,
                outgoing,
            )

            # ignoring nocaustics for now

            # handle opacity
            # if !isapprox(materialPoint.opacity, 1) &&
            if materialPoint.opacity < 1 &&
               rand(Float32) >= materialPoint.opacity
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

            # set hit variables
            # if (bounce == 0)
            #     hit = true
            # end

            # accumulate emission
            radiance +=
                weight *
                (dot(normal, outgoing) >= 0 ? materialPoint.emission : zeroSV3f)

            # next direction
            incoming = zeroSV3f
            # if materialPoint.roughness != 0
            if !isDelta(materialPoint)
                # BSDF materials
                if rand(Float32) < 0.5f0
                    incoming = sampleBSDF(materialPoint, normal, outgoing)
                else
                    incoming = sampleLights(scene, lights, position)
                end

                if incoming == zeroSV3f
                    break
                end
                # weight .*
                weight *=
                    evalBSDF(materialPoint, normal, outgoing, incoming) / (
                        0.5f0 *
                        pdfBSDF(materialPoint, normal, outgoing, incoming) +
                        0.5f0 *
                        pdfLights(scene, bvh, lights, position, incoming)
                    )
            else
                # Delta materials
                incoming = sampleDelta(materialPoint, normal, outgoing)
                # TODO: understand if this if is needed
                # if incoming == zeroSV3f
                #     break
                # end
                weight *=
                    evalDelta(materialPoint, normal, outgoing, incoming) /
                    pdfDelta(materialPoint, normal, outgoing, incoming)
            end

            # update volume stack
            if (
                isVolumetric(materialPoint) &&
                (dot(normal, outgoing) * dot(normal, incoming)) < 0
            )
                if (isempty(volumeStack))
                    materialPoint = evalMaterial(scene, instance, intersection)
                    push!(volumeStack, materialPoint)
                else
                    pop!(volumeStack)
                end
            end

            # setup next iteration
            newRay = Ray(position, incoming)

        else
            outgoing = -newRay.direction
            position = ray.origin + ray.direction * intersection.distance
            vsdf = volumeStack[end]

            incoming = zeroSV3f
            if rand(Float32) < 0.5f0
                incoming = sampleScattering(vsdf, outgoing)
            else
                incoming = sampleLights(scene, lights, position)
            end

            if incoming == zeroSV3f
                break
            end

            weight *=
                evalScattering(vsdf, outgoing, incoming) / (
                    0.5f0 * sampleScatteringPdf(vsdf, outgoing, incoming) +
                    0.5f0 * pdfLights(scene, bvh, lights, position, incoming)
                )

            # setup next iteration
            newRay = Ray(position, incoming)
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
    end
    return radiance
end

@inline function isVolumetric(material::MaterialPoint)
    return material.type == "refractive" ||
           material.type == "volume" ||
           material.type == "subsurface"
end

@inline function isDelta(material::MaterialPoint)
    return (material.type == "reflective" && material.roughness == 0) ||
           (material.type == "refractive" && material.roughness == 0) ||
           (material.type == "transparent" && material.roughness == 0) ||
           (material.type == "volume")
end