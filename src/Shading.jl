module Shading

using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot

export shaderEyelight,
    shaderMaterial, shaderNormal, shaderIndirectNaive, shaderEyelightBsdf

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
    normal = evalNormal(shape, intersection, frame, outgoing)

    radiance::SVec3f = normal * 0.5 .+ 0.5

    return radiance
end

function shaderEyelightBsdf(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f

    # initialize
    radiance::SVec3f = (0, 0, 0)
    weight::SVec3f = (1, 1, 1)
    maxBounces = 4
    opacityBounce = 0

    for bounce = 1:maxBounces

        # intersect 
        intersection::Intersection = intersectScene(ray, scene, bvh, false)
        if !intersection.hit
            radiance += weight .* evalEnvironment(scene, ray.direction)
            break
        end

        # extract intersection data
        instance::Instance = scene.instances[intersection.instanceIndex]
        shape::Shape = scene.shapes[instance.shapeIndex]
        frame::Frame = instance.frame

        # eval position
        outgoing::SVec3f = -ray.direction
        positon::SVec3f = evalPosition(
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
        shapeColor::SVec4f = (1, 1, 1, 1)

        # evaluate color and emission of the  material
        material::Material = scene.materials[instance.materialIndex]

        # emission
        materialEmissionTex::SVec4f =
            evalTexture(scene, material.emissionTex, textureX, textureY, true)
        materialEmission::SVec3f =
            material.emission .* xyz(materialEmissionTex) .* xyz(shapeColor)

        # color
        materialColorTex::SVec4f =
            evalTexture(scene, material.colorTex, textureX, textureY, true)
        materialColor::SVec3f =
            material.color .* xyz(materialColorTex) .* xyz(shapeColor)

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
        radiance += weight .* emission

        # brdf + light
        radiance +=
            (weight * pi) .* evalBsdfCos(
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

function evalBsdfCos(
    materialType::String,
    materialRoughness::Float32,
    materialColor::SVec3f,
    normal,
    outgoing,
    incoming,
)::SVec3f
    if materialRoughness == 0
        return SVec3f(0, 0, 0)
        error("wait what")
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
        material.emission .*
        xyz(evalTexture(scene, material.emissionTex, textureX, textureY))
    color::SVec3f =
        material.color .*
        xyz(evalTexture(scene, material.colorTex, textureX, textureY))

    # materialColor = evalMaterialColor(scene, intersection)

    # radiance::SVec3f = emission + abs(dot(normal, outgoing)) .* color
    radiance::SVec3f = abs(dot(normal, outgoing)) .* color

    return SVec3f(radiance.x, radiance.y, radiance.z)
end

function shaderIndirectNaive(
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    bounce::Int = 1,
)::SVec3f
    maxBounce = 10

    # INTERSECT SCENE
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

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
    normal = evalNormal(shape, intersection, frame)
    textureX, textureY = evalTexcoord(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )
    outgoing = -ray.direction

    # NORMAL CORRECTIONS
    if !isempty(shape.triangles) && dot(normal, outgoing) < 0
        normal = -normal
    end

    # EVALUATE MATERIAL
    material::Material = scene.materials[instance.materialIndex]
    materialEmissionTex::SVec4f =
        evalTexture(scene, material.emissionTex, textureX, textureY, true)
    materialColorTex::SVec4f =
        evalTexture(scene, material.colorTex, textureX, textureY, true)
    emission::SVec3f = material.emission .* xyz(materialEmissionTex)
    color::SVec3f = material.color .* xyz(materialColorTex)
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
    radiance::SVec3f = emission

    # EXIT IF RAY IS DONE
    if bounce >= maxBounce
        return radiance
    end

    # COMPUTE ILLUMINATION

    # everything matte for now
    incoming::SVec3f = sampleHemisphereCos(normal)
    lighting::SVec3f =
        shaderIndirectNaive(scene, Ray(position, incoming), bvh, bounce + 1)
    radiance += color .* lighting

    # if incoming == SVec3f(0, 0, 0)
    #     return radiance
    # end
    # if dot(normal, incoming) * dot(normal, outgoing) > 0
    #     radiance +=
    #         materialColor .*
    #         shaderIndirectNaive(scene, Ray(position, incoming), bvh, bounce + 1)
    # end
    # return radiance
    return radiance
end

function shaderMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    radiance = SVec3f(0, 0, 0)
    maxBounce = 2
    newRay = ray

    for b = 1:maxBounce
        weight::Float32 = 1.0f0 / b
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance = evalEnvironment(scene, ray.direction)
            return weight * radiance
        end

        # compute normal of the point hit
        instance::Instance = scene.instances[intersection.instanceIndex]
        frame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]
        material::Material = scene.materials[instance.materialIndex]

        outgoing = -ray.direction
        position = evalShadingPosition(scene, intersection, outgoing)
        normal = evalNormal(shape, intersection, frame, outgoing)
        materialColor = evalMaterialColor(scene, intersection)

        radiance = evalEmission(material, normal, outgoing) .* materialColor

        if material.type == "matte"
            incoming = sampleHemisphereCos(normal)
            if incoming == SVec3f(0, 0, 0)
                break
            end
            if dot(normal, incoming) * dot(normal, outgoing) > 0
                radiance = linInterp(
                    radiance,
                    radiance / pi .* abs(dot(normal, outgoing)),
                    weight,
                )
            end
        elseif material.type == "reflective"
            if material.roughness == 0
                incoming = reflect(outgoing, normal)
                if dot(normal, incoming) * dot(normal, outgoing) <= 0
                    return weight * radiance
                end
                radiance = linInterp(
                    radiance,
                    radiance .* fresnelSchlick(materialColor, normal, incoming),
                    weight,
                )
            else
                exponent = 2.0f0 / material.roughness^2
                halfway = sampleHemisphereCosPower(exponent, normal)
                incoming = reflect(outgoing, halfway)
                if dot(normal, incoming) * dot(normal, outgoing) <= 0
                    return weight * radiance
                end
                up_normal = dot(normal, outgoing) <= 0 ? -normal : normal
                halfway = norm(incoming + outgoing)
                F = fresnelConductor(
                    reflectivityToEta(materialColor),
                    SVec3f(0, 0, 0),
                    halfway,
                    incoming,
                )
                D = microfacetDistribution(
                    material.roughness,
                    up_normal,
                    halfway,
                )
                G = microfacetShadowing(
                    material.roughness,
                    up_normal,
                    halfway,
                    outgoing,
                    incoming,
                )
                radiance = linInterp(
                    radiance,
                    radiance .* F .* D .* G ./ (
                        4 .* dot(up_normal, outgoing) .*
                        dot(up_normal, incoming)
                    ) .* abs(dot(up_normal, incoming)),
                    weight,
                )
            end
            # elseif material.roughness == 0

        else
            incoming = sampleHemisphereCos(normal)
            if incoming == SVec3f(0, 0, 0)
                break
            end
            #if dot(normal, incoming) * dot(normal, outgoing) > 0
            radiance = linInterp(
                radiance,
                materialColor / pi .* abs(dot(normal, incoming)),
                weight,
            )
            #end
        end

        newRay = Ray(position, incoming)

        # radiance = 0.5 .* (normal .+ 1) .* color
        #radiance::SVec3f = abs(dot(normal, outgoing)) .* radiance
    end
    return radiance
end

@inline function evalEmission(
    material::Material,
    normal::SVec3f,
    outgoing::SVec3f,
)
    ifelse(dot(normal, outgoing) >= 0, material.emission, SVec3f(0, 0, 0))
end

@inline function fresnelSchlick(
    specular::SVec3f,
    normal::SVec3f,
    outgoing::SVec3f,
)::SVec3f
    if specular == SVec3f(0, 0, 0)
        return SVec3f(0, 0, 0)
    end
    cosine = dot(normal, outgoing)
    specular .+
    (1.0f0 .- specular) .* clamp(1.0f0 .- abs(cosine), 0.0f0, 1.0f0)^5.0f0
end

@inline function reflectivityToEta(color::SVec3f)::SVec3f
    reflectivity::SVec3f = clamp.(color, 0.0f0, 1.0f0)
    return (1.0f0 .+ sqrt.(reflectivity)) ./ (1.0f0 .- sqrt.(reflectivity))
end

@inline function reflect(w::SVec3f, n::SVec3f)::SVec3f
    return -w + 2.0f0 * dot(w, n) * n
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
        # maybe this is wrong because of the + 1 with no parentesis?
        return roughness2 / ((pi * (cosine2 * roughness2 + 1 - cosine2)^2))
    else
        return exp((cosine2 - 1) / (roughness2 * cosine2)) /
               (pi * roughness2 * cosine2 * cosine2)
    end
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

function sampleHemisphere(normal::SVec3f)::SVec3f
    ruvx = rand()
    z = rand()  # ruvy
    r = sqrt(clamp(1 - z * z, 0.0f0, 1.0f0))
    phi = 2 * pi * ruvx
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

function sampleHemisphereCos(normal::SVec3f)::SVec3f
    ruvx = rand()
    z = sqrt(rand())  # ruvy
    r = sqrt(1 - z * z)
    phi = 2 * pi * ruvx
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
end

@inline function sampleHemisphereCosPower(
    exponent::Float32,
    normal::SVec3f,
)::SVec3f
    ruxy = rand()
    ruvy = rand()
    z = ruvy^(1.0f0 / (exponent + 1.0f0))
    r = sqrt(1.0f0 - z * z)
    phi = 2.0f0 * pi * ruxy
    localDirection = SVec3f(r * cos(phi), r * sin(phi), z)
    return transformDirection(basisFromz(normal), localDirection)
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