
using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot, cross

export shadeEyelight, shadeNormal, shadeColor

# shader using the color of the material as false color
function shadeColor(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    instance::Instance = scene.instances[intersection.instanceIndex]
    material::Material = scene.materials[instance.materialIndex]

    radiance = material.color

    return radiance
end

# shader using normal orientation as false color
function shadeNormal(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh)

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

# Eyelight shader
function shadeEyelight(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f

    # initialize
    radiance::SVec3f = zeroSV3f
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
        emission = dot(normal, outgoing) >= 0 ? materialEmission : zeroSV3f
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

# function shaderEyelight(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
#     intersection::Intersection = intersectScene(ray, scene, bvh, false)

#     if !intersection.hit
#         radiance = evalEnvironment(scene, ray.direction)
#         return radiance
#     end

#     # EVALUATE GEOMETRY 
#     instance::Instance = scene.instances[intersection.instanceIndex]
#     frame::Frame = instance.frame
#     shape::Shape = scene.shapes[instance.shapeIndex]
#     normal = evalNormal(shape, intersection, frame)

#     textureX, textureY = evalTexcoord(
#         scene,
#         instance,
#         intersection.elementIndex,
#         intersection.u,
#         intersection.v,
#     )
#     outgoing = -ray.direction

#     # EVALUATE MATERIAL
#     material::Material = scene.materials[instance.materialIndex]
#     emission::SVec3f =
#         material.emission *
#         xyz(evalTexture(scene, material.emissionTex, textureX, textureY))
#     color::SVec3f =
#         material.color *
#         xyz(evalTexture(scene, material.colorTex, textureX, textureY))

#     # materialColor = evalMaterialColor(scene, intersection)

#     # radiance::SVec3f = emission + abs(dot(normal, outgoing)) * color
#     radiance::SVec3f = abs(dot(normal, outgoing)) * color

#     return SVec3f(radiance.x, radiance.y, radiance.z)
# end
