module Shading

using ..World
using ..Intersect
using ..Eval
using ..Bvh
using ..Types
using ..Algebra

using StaticArrays: dot

export shaderEyelight,
    shaderMaterial, shaderNormal, shaderIndirectNaive, shaderIndirect

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

function shaderEyelight(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]
    material::Material = scene.materials[instance.materialIndex]

    outgoing = -ray.direction
    normal = evalNormal(shape, intersection, frame, outgoing)

    materialColor = evalMaterialColor(scene, intersection)

    radiance::SVec3f = abs(dot(normal, outgoing)) .* materialColor

    return radiance
end

function shaderIndirectNaive(
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    bounce::Int = 1,
)::SVec3f
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    maxBounce = 6
    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]
    material::Material = scene.materials[instance.materialIndex]

    outgoing = -ray.direction
    position = evalShadingPosition(scene, intersection, outgoing)
    normal = evalNormal(shape, intersection, frame, outgoing)
    materialColor = evalMaterialColor(scene, intersection)
    radiance = material.emission
    if bounce >= maxBounce
        return radiance
    end

    incoming = sampleHemisphereCos(normal)
    if incoming == SVec3f(0, 0, 0)
        return radiance
    end
    if dot(normal, incoming) * dot(normal, outgoing) > 0
        radiance +=
            materialColor .*
            shaderIndirectNaive(scene, Ray(position, incoming), bvh, bounce + 1)
    end
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