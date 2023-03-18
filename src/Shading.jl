module Shading

using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot, cross

export shadePath, shadeMaterial

include("Baseshaders.jl")
include("Pathtrace.jl")
include("Raytrace.jl")

# this is meant to be a copy of trace_naive but with only matte
function shadePath(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    radiance::SVec3f = zeroSV3f
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
            radiance =
                muladd.(weight, evalEnvironment(scene, ray.direction), radiance)
            #radiance += weight * evalEnvironment(scene, ray.direction)
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
        radiance =
            muladd.(
                weight,
                evalEmission(materialEmission, normal, outgoing),
                radiance,
            )
        #radiance += weight * evalEmission(materialEmission, normal, outgoing)

        # next direction
        incoming::SVec3f = zeroSV3f
        if materialRoughness != 0
            incoming = sampleBsdfCos(
                material.type,
                materialRoughness,
                materialColor,
                normal,
                outgoing,
            )
            if incoming == zeroSV3f
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
        if weight == zeroSV3f || !all(isfinite.(weight))
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
        return zeroSV3f
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
        return zeroSV3f
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
        return zeroSV3f
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
        return zeroSV3f
    end
    return materialColor / pi * abs(dot(normal, incoming))
end

function evalEmission(materialEmission, normal, outgoing)::SVec3f
    return dot(normal, outgoing) >= 0 ? materialEmission : zeroSV3f
end

end