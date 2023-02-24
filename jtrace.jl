using Images
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using .Types
using .Algebra
using .Bvh
using .Intersect
using .World
using .Eval

# main entry point to the program
function run(
    scenePath::String,
    shader::Function = shaderEyelight,
    width = 1920,
    numSamples = 2,
)
    # generate scene
    scene = loadJsonScene(joinpath(baseDir, scenePath))
    # return scene

    camera = scene.cameras[1]
    height = Int(round(width / camera.aspect))

    # build bvh
    bvh = makeSceneBvh(scene)

    # generate empty starting image
    image = zeros(SVec3f, height, width)

    # call the function to trace samples
    traceSamples(shader, image, scene, width, height, numSamples, bvh, camera)

    # save the resulting image
    rgbImage = zeros(RGB, size(image))
    Threads.@threads for i = 1:size(image)[1]
        Threads.@threads for j = 1:size(image)[2]
            rgbImage[i, j] = RGB(image[i, j]...)
        end
    end
    save("out/prova.png", rgbImage)
end

function traceSamples(
    shader,
    image,
    scene,
    imwidth,
    imheight,
    numSamples,
    bvh,
    camera,
)
    for s = 1:numSamples
        Threads.@threads for i = 1:size(image)[2]
            Threads.@threads for j = 1:size(image)[1]
                color = traceSample(
                    shader,
                    i,
                    j,
                    scene,
                    camera,
                    imwidth,
                    imheight,
                    bvh,
                )

                weight::Float32 = 1 / s
                image[j, i] =
                    clamp.(linInterp(image[j, i], color, weight), 0.0f0, 1.0f0)
            end
        end
    end
end

function traceSample(
    shader::Function,
    i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int,
    bvh,
)::SVec3f

    # send a ray
    ray = sampleCamera(camera, i, j, imwidth, imheight)

    # call the shader
    radiance = shader(scene, ray, bvh)

    return radiance
end

function shaderColor(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    intersection::Intersection = intersectScene(ray, scene)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    instance::Instance = scene.instances[intersection.instanceIndex]
    material::Material = scene.materials[instance.materialIndex]

    radiance = material.color

    return radiance
end

function shaderNormal(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]

    normal = evalNormal(shape, intersection, frame)

    radiance::SVec3f = normal * 0.5 .+ 0.5

    return radiance
end

function shaderEyelight(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    intersection::Intersection = intersectScene(ray, scene, bvh, false)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        # radiance = SVec3f(0, 0, 0)
        return radiance
    end

    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]
    material::Material = scene.materials[instance.materialIndex]

    normal = evalNormal(shape, intersection, frame)

    outgoing = -ray.direction

    materialColor = evalMaterialColor(scene, intersection)

    radiance::SVec3f = abs(dot(normal, outgoing)) .* materialColor

    return radiance
end

function evalMaterialColor(scene::Scene, intersection::Intersection)::SVec3f
    instance::Instance = scene.instances[intersection.instanceIndex]
    material::Material = scene.materials[instance.materialIndex]
    elementIndex::Int = intersection.elementIndex
    u::Float32 = intersection.u
    v::Float32 = intersection.v

    textureU, textureV = evalTexcoord(scene, instance, elementIndex, u, v)

    colorTexture = evalTexture(scene, material.colorTex, textureU, textureV)

    pointColor = material.color .* xyz(colorTexture)

    return pointColor
end

function evalTexcoord(
    scene::Scene,
    instance::Instance,
    elementIndex::Int,
    u::Float32,
    v::Float32,
)
    shape::Shape = scene.shapes[instance.shapeIndex]
    if isempty(shape.textureCoords)
        return u, v
    end
    if !isempty(shape.triangles)
        t1, t2, t3 = @view shape.triangles[elementIndex, :]
        return interpolateTriangle(
            SVec2f(@view shape.textureCoords[t1, :]),
            SVec2f(@view shape.textureCoords[t2, :]),
            SVec2f(@view shape.textureCoords[t3, :]),
            u,
            v,
        )

    elseif !isempty(shape.quads)
        q1, q2, q3, q4 = @view shape.quads[elementIndex, :]
        return interpolateQuad(
            SVec2f(@view shape.textureCoords[q1, :]),
            SVec2f(@view shape.textureCoords[q2, :]),
            SVec2f(@view shape.textureCoords[q3, :]),
            SVec2f(@view shape.textureCoords[q4, :]),
            u,
            v,
        )
    else
        error("No triangles or quads in this shape")
    end
end

function evalShadingPosition(
    scene::Scene,
    intersection::Intersection,
    outgoing::SVec3f,
)::SVec3f
    instance = scene.instances[intersection.instanceIndex]
    element = intersection.elementIndex
    u = intersection.u
    v = intersection.v
    shape = scene.shapes[instance.shapeIndex]
    if (!isempty(shape.triangles) || !isempty(shape.quads))
        return evalPosition(scene, instance, element, u, v)
    end
end

function evalPosition(
    scene::Scene,
    instance::Instance,
    element::Int,
    u::Float32,
    v::Float32,
)
    shape = scene.shapes[instance.shapeIndex]
    if !isempty(shape.triangles)
        t1, t2, t3 = @view shape.triangles[element, :]
        return transformPoint(
            instance.frame,
            interpolateTriangle(
                SVec3f(@view shape.positions[t1, :]),
                SVec3f(@view shape.positions[t2, :]),
                SVec3f(@view shape.positions[t3, :]),
                u,
                v,
            ),
        )
    elseif !isempty(shape.quads)
        q1, q2, q3, q4 = @view shape.quads[element, :]
        return transformPoint(
            instance.frame,
            interpolateQuad(
                SVec3f(@view shape.positions[q1, :]),
                SVec3f(@view shape.positions[q2, :]),
                SVec3f(@view shape.positions[q3, :]),
                SVec3f(@view shape.positions[q4, :]),
                u,
                v,
            ),
        )
    else
        error("No triangles or quads in this shape")
    end
end

function shaderMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    radiance = SVec3f(0, 0, 0)
    maxBounce = 2
    newRay = ray

    for b = 1:maxBounce
        intersection::Intersection = intersectScene(newRay, scene, bvh, false)

        if !intersection.hit
            radiance = evalEnvironment(scene, newRay.direction)
            return radiance
        end

        # compute normal of the point hit
        instance::Instance = scene.instances[intersection.instanceIndex]
        frame::Frame = instance.frame
        shape::Shape = scene.shapes[instance.shapeIndex]
        material::Material = scene.materials[instance.materialIndex]

        outgoing = -newRay.direction
        position = evalShadingPosition(scene, intersection, outgoing)
        normal = evalNormal(shape, intersection, frame)
        materialColor = evalMaterialColor(scene, intersection)

        radiance = material.emission

        if material.type == "matte"
            incoming = sampleHemisphereCos(normal)
            if dot(normal, incoming) * dot(normal, outgoing) > 0
                radiance += materialColor / pi .* abs(dot(normal, outgoing))
            end
        elseif material.type == "reflective"
            if material.roughness == 0
                incoming = reflect(outgoing, normal)
                if dot(normal, incoming) * dot(normal, outgoing) <= 0
                    return radiance
                end
                radiance +=
                    materialColor .*
                    fresnelSchlick(materialColor, normal, incoming)
            else
                exponent = 2 / material.roughness^2
                halfway = sampleHemisphereCosPower(exponent, normal)
                incoming = reflect(outgoing, halfway)
                if dot(normal, incoming) * dot(normal, outgoing) <= 0
                    return radiance
                end
                up_normal = dot(normal, outgoing) <= 0 ? -normal : normal
                halfway = normalize(incoming + outgoing)
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
                radiance +=
                    materialColor .* F .* D .* G ./ (
                        4 .* dot(up_normal, outgoing) .*
                        dot(up_normal, incoming)
                    ) .* abs(dot(up_normal, incoming))
            end
            # elseif material.roughness == 0

        else
            incoming = sampleHemisphereCos(normal)
            if dot(normal, incoming) * dot(normal, outgoing) > 0
                radiance += materialColor / pi .* abs(dot(normal, incoming))
            end
        end

        newRay = Ray(position, incoming)

        # radiance = 0.5 .* (normal .+ 1) .* color
        #radiance::SVec3f = abs(dot(normal, outgoing)) .* radiance
    end
    return radiance
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
    r = sqrt(1 - z * z) # maybe clamp here?
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
    z = normalize(v)
    sign = copysign(1.0f0, z.z)
    a = -1.0f0 / (sign + z.z)
    b = z.x * z.y * a
    x = SVec3f(1.0f0 + sign * z.x * z.x * a, sign * b, -sign * z.x)
    y = SVec3f(b, sign + z.y * z.y * a, -z.y)
    return Frame(x, y, z)
end

function evalEnvironment(scene::Scene, direction::SVec3f)::SVec3f
    # background = SVec3f(0.105, 0.443, 0.90)
    emission = SVec3f(0, 0, 0)
    for env in scene.environments
        emission += evalEnvironment(scene, env, direction)
    end
    return emission
end

function evalEnvironment(
    scene::Scene,
    env::Environment,
    direction::SVec3f,
)::SVec3f
    wl::SVec3f = transformDirection(inverse(env.frame), direction)
    textureX::Float32 = atan(wl[3], wl[1]) / (2 * pi)
    textureY::Float32 = acos(clamp(wl[2], -1, 1)) / pi

    if textureX < 0
        textureX += 1
    end

    return env.emission .*
           xyz(evalTexture(scene, env.emissionTex, textureX, textureY))
end

function evalTexture(
    scene::Scene,
    textureIdx::Int,
    textureX::Float32,
    textureY::Float32,
)::SVec4f
    if textureIdx == -1 || textureIdx == 0
        return SVec4f(1, 1, 1, 1)
    end

    texture = scene.textures[textureIdx]
    return evalTexture(texture, textureX, textureY)
end

function evalTexture(
    texture::Texture,
    textureX::Float32,
    textureY::Float32,
)::SVec4f
    if !isempty(texture.image)
        sizeX, sizeY = size(texture.image)
    elseif !isempty(texture.hdrImage)
        (texture.hdrImage)
        sizeX, sizeY = size(texture.hdrImage)
    else
        error("Texture contains no image")
    end

    asLinear = false
    clampToEdge = texture.clamp
    # noInterpolation = texture.nearest
    noInterpolation = true
    s = 0.0f0
    t = 0.0f0

    if clampToEdge
        s = clamp(textureX, 0, 1.0f0) * sizeX
        t = clamp(textureY, 0, 1.0f0) * sizeY
    else
        s = rem(textureX, 1.0f0) * sizeX
        if (s <= 0)
            s += sizeX
        end
        t = rem(textureY, 1.0f0) * sizeY
        if (t <= 0)
            t += sizeY
        end
    end

    i::Int = clamp(Int(ceil(s)), 1, sizeX)
    j::Int = clamp(Int(ceil(t)), 1, sizeY)

    ii::Int = (i + 1) % sizeX
    jj::Int = (j + 1) % sizeY
    u::Float32 = s - i
    v::Float32 = t - j

    if noInterpolation
        return lookupTexture(texture, i, j)
    else
        return (
            lookupTexture(texture, i, j) * (1 - u) * (1 - v) +
            lookupTexture(texture, i, jj) * (1 - u) * v +
            lookupTexture(texture, ii, j) * u * (1 - v) +
            lookupTexture(texture, ii, jj) * u * v
        )
    end
end

function lookupTexture(texture::Texture, i::Int, j::Int)::SVec4f
    if !isempty(texture.image)
        rgba = texture.image[i, j]
        color = SVec4f(rgba.r, rgba.g, rgba.b, rgba.alpha)
    elseif !isempty(texture.hdrImage)
        rgb = texture.hdrImage[i, j]
        color = SVec4f(rgb.r, rgb.g, rgb.b, 1)
    else
        error("Texture contains no image")
    end
    return color
end

function evalNormalSphere(ray::Ray, sphereCenter::SVec3f)
    # compute coordinates of point hit
    pointHit::SVec3f = ray.origin + ray.tmin * ray.direction
    # compute normal: n = (p-c) / |p-c|
    normal = unitVector(pointHit - sphereCenter)
    return normal
end
