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
    radiance = SVec3f(0, 0, 0)
    maxBounce = 128
    newRay = ray
    opbounce::Int8 = 0
    weight = SVec3f(1, 1, 1)

    for bounce = 1:maxBounce
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
            weight * (
                dot(normal, outgoing) >= 0 ? materialPoint.emission :
                SVec3f(0, 0, 0)
            )

        incoming = SVec3f(0, 0, 0)

        # next direction
        if materialPoint.roughness != 0
            # BSDF materials
            incoming = sampleBSDF(materialPoint, normal, outgoing)
            if incoming == SVec3f(0, 0, 0)
                break
            end
            weight *=
                evalBSDF(materialPoint, normal, outgoing, incoming) /
                pdfBSDF(materialPoint, normal, outgoing, incoming)
        else
            # Delta materials
            incoming = sampleDelta(materialPoint, normal, outgoing)
            if incoming == SVec3f(0, 0, 0)
                break
            end
            weight *=
                evalDelta(materialPoint, normal, outgoing, incoming) /
                pdfDelta(materialPoint, normal, outgoing, incoming)
        end

        # check weight
        if (weight == SVec3f(0, 0, 0) || !all(isfinite.(weight)))
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

# TODO: move to Eval.jl
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