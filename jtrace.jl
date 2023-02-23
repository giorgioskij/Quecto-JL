using Images
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using .Types
using .Algebra
using .Bvh

const global shapeBvhDepth = 25
const global masterBvhDepth = 19
#const global shapeBvhDepth = 20
#const global masterBvhDepth = 3

# # (ceil(log2(count((n) -> !n.internal, bvh.nodes))) + 2) * 2
# const global masterBvhDepth =
#     ceil(log2(size(bvh.bvh.primitives, 1)) + 2) * 2
# # ( ( [count((n) -> !n.internal, s.bvh.nodes) for s in bvh.shapes] |> max |> log2) + 1) |> ceil * 2
# const global shapeBvhDepth =
#     ceil(log2(max([size(s.bvh.primitives, 1) for s in bvh.shapes])) + 2) *

# main entry point to the program
function run(
    scenePath::String,
    shader::Function = shaderMaterial,
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
    # camera = scene.cameras[1]
    # loop over pixels
    # println("Starting creation of image...")
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

    color = material.color

    # radiance = 0.5 .* (normal .+ 1) .* color
    radiance::SVec3f = abs(dot(normal, outgoing)) .* material.color

    return radiance
end

function shaderMaterial(scene::Scene, ray::Ray, bvh::SceneBvh)::SVec3f
    # intersection::Intersection = intersectScene(ray, scene)
    radiance = SVec3f(0, 0, 0)
    maxBounce = 4

    for b = 1:maxBounce
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

        radiance = material.emission
        normal = evalNormal(shape, intersection, frame)

        # do positions
        #position = evalPosition(shape, intersection, frame)
        outgoing = -ray.direction

        if material.type == "matte"
            incoming = sampleHemisphereCos(normal)
            if dot(normal, incoming) * dot(normal, outgoing) <= 0
                return SVec3f(0, 0, 0)
            end
            radiance += material.color / pi .* abs(dot(normal, outgoing))
        elseif material.type == "reflective"
            if material.roughness == 0
                incoming = reflect(outgoing, normal)
                # if dot(normal, incoming) * dot(normal, outgoing) <= 0
                #     return SVec3f(0, 0, 0)
                # end
                radiance +=
                    material.color .*
                    fresnelSchlick(material.color, normal, incoming)
            else
                exponent = 2 / material.roughness^2
                halfway = sampleHemisphereCosPower(exponent, normal)
                incoming = reflect(outgoing, halfway)
                if dot(normal, incoming) * dot(normal, outgoing) <= 0
                    return SVec3f(0, 0, 0)
                end
                up_normal = dot(normal, outgoing) <= 0 ? -normal : normal
                halfway = normalize(incoming + outgoing)
                F = fresnelConductor(
                    reflectivityToEta(material.color),
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
                    material.color .* F .* D .* G ./ (
                        4 .* dot(up_normal, outgoing) .*
                        dot(up_normal, incoming)
                    ) .* abs(dot(up_normal, incoming))
            end

            # elseif material.roughness == 0

        else
            radiance = abs(dot(normal, outgoing)) .* material.color
        end

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
    return w + 2.0f0 * dot(w, n) * n
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

function basisFromz(v::SVec3f)
    z = normalize(v)
    sign = copysign(1.0f0, z.z)
    a = -1.0f0 / (sign + z.z)
    b = z.x * z.y * a
    x = SVec3f(1.0f0 + sign * z.x * z.x * a, sign * b, -sign * z.x)
    y = SVec3f(b, sign + z.y * z.y * a, -z.y)
    return Frame(x, y, z)
end

function intersectTriangle(
    ray::Ray,
    triangle::Triangle,
    instanceIndex::Int64,
    triangleIndex::Int64,
)::Intersection
    edge1 = triangle.y - triangle.x
    edge2 = triangle.z - triangle.x
    pvec = cross(ray.direction, edge2)
    det = dot(edge1, pvec)

    if det == 0
        return Intersection(false)
    end

    inverseDet = 1.0f0 / det

    tvec = ray.origin - triangle.x
    u = dot(tvec, pvec) * inverseDet
    if u < 0 || u > 1
        return Intersection(false)
    end

    qvec = cross(tvec, edge1)
    v = dot(ray.direction, qvec) * inverseDet
    if v < 0 || u + v > 1
        return Intersection(false)
    end

    t = dot(edge2, qvec) * inverseDet
    if t < ray.tmin || t > ray.tmax
        return Intersection(false)
    end

    return Intersection(true, instanceIndex, triangleIndex, u, v, t)
end

function intersectQuad(
    ray::Ray,
    quad::Quad,
    instanceIndex::Int64,
    quadIndex::Int64,
)
    if (quad.c == quad.d)
        return intersectTriangle(
            ray,
            Triangle(quad.a, quad.b, quad.c),
            instanceIndex,
            quadIndex,
        )
    end

    isec1 = intersectTriangle(
        ray,
        Triangle(quad.a, quad.b, quad.d),
        instanceIndex,
        quadIndex,
    )
    if (isec1.hit)
        isec1 = Intersection(
            true,
            instanceIndex,
            quadIndex,
            isec1.u,
            isec1.v,
            isec1.distance,
        )
    end
    isec2 = intersectTriangle(
        ray,
        Triangle(quad.c, quad.d, quad.b),
        instanceIndex,
        quadIndex,
    )
    if (isec2.hit)
        isec2 = Intersection(
            true,
            instanceIndex,
            quadIndex,
            1 - isec2.u,
            1 - isec2.v,
            isec2.distance,
        )
    end

    if (isec1.hit && !isec2.hit)
        return isec1
    elseif (isec2.hit && !isec1.hit)
        return isec2
    elseif (isec1.hit && isec2.hit)
        return isec1.distance < isec2.distance ? isec1 : isec2
    else
        return isec1
    end
end

function intersectScene(
    ray::Ray,
    scene::Scene,
    sceneBvh::SceneBvh,
    findAny::Bool,
)::Intersection
    masterBvh = sceneBvh.bvh

    if (isempty(masterBvh.nodes))
        return Intersection(false)
    end

    # node stack
    nodeCur = 1
    # PROF: this is slow - 3200 ms in @profview (1920, 1080, 2)
    # nodeStack = zeros(Int, 128)
    # let's try with mvector
    # with MVector the creation is much faster, at 1800 ms 
    # With Int16 the speedup is substantial 4.5 --> 5.1 seconds
    # but UInt16 may not be enough, UInt32 should be
    # nodeStack = zeros(MVector{128,UInt32})
    nodeStack = MVector{masterBvhDepth,UInt32}(undef)

    # this is very esoteric optimization, be careful with side effects
    # lets pre allocate this guys to pass to the calls of intersectShapeBvh
    preallocatedNodeStackForMySon = MVector{shapeBvhDepth,UInt32}(undef)

    nodeStack[nodeCur] = 1
    nodeCur += 1

    # init intersection
    intersection::Intersection = Intersection(false)

    rayDInv = SVec3f(
        1.0f0 / ray.direction.x,
        1.0f0 / ray.direction.y,
        1.0f0 / ray.direction.z,
    )
    rayDSign = SVec3i(rayDInv.x < 0, rayDInv.y < 0, rayDInv.z < 0)

    # walking stack
    while (nodeCur != 1)
        # if nodeCur > maxNodeCur
        #     println("nodecur")
        # end

        # grab node
        nodeCur -= 1
        node = masterBvh.nodes[nodeStack[nodeCur]]

        # intersect bbox
        if !intersectBbox(ray, rayDInv, node.bbox)
            continue
        end

        # intersect node, switching based on node type
        # for each type, iterate over the primitive list
        if node.internal
            # for internal nodes, attempts to proceed along the
            # split axis from smallest to largest nodes
            if (rayDSign[node.axis] != 0)
                nodeStack[nodeCur] = node.start + 0
                nodeCur += 1
                nodeStack[nodeCur] = node.start + 1
                nodeCur += 1
            else
                nodeStack[nodeCur] = node.start + 1
                nodeCur += 1
                nodeStack[nodeCur] = node.start + 0
                nodeCur += 1
            end
        else
            for idx = node.start:node.start+node.num-1
                instance = scene.instances[masterBvh.primitives[idx]]
                invRay = transformRay(inverse(instance.frame, true), ray)

                # to understand this part because we create a lot of objects and then 
                # we possibly rewrite them, we need to understand the performance impact
                # and maybe create the objects only at the end of the loop?
                sIntersection = intersectShapeBvh!(
                    sceneBvh.shapes[instance.shapeIndex],
                    scene.shapes[instance.shapeIndex],
                    invRay,
                    findAny,
                    preallocatedNodeStackForMySon,
                )
                if (!sIntersection.hit)
                    continue
                end
                intersection = Intersection(
                    true,
                    masterBvh.primitives[idx],
                    sIntersection.elementIndex,
                    sIntersection.u,
                    sIntersection.v,
                    sIntersection.distance,
                )
                # change ray tmax
                ray = Ray(
                    ray.origin,
                    ray.direction,
                    ray.tmin,
                    sIntersection.distance,
                )
            end
        end
        if (findAny && intersection.hit)
            return intersection
        end
    end
    return intersection
end

function intersectShapeBvh!(
    shapeBvh::ShapeBvh,
    shape::Shape,
    ray::Ray,
    findAny::Bool,
    nodeStack::MVector, # this is passed to avoid multiple initializations
)::ShapeIntersection
    bvh = shapeBvh.bvh

    # check empty
    if isempty(bvh.nodes)
        return ShapeIntersection(false)
    end

    # node stack
    # nodeStack = zeros(Int, 128)
    # PROF: same thing here. we need MVectors
    # nodeStack = Vector{Int}(undef, 128)
    # just going from vector to Mvector we speed up, but smaller ints are even better
    # UInt16 is definitely not enough :(, we get truncation error
    # OMFG GOING FROM Int64 to UInt32 we go from 4.05 seconds to 2.8!!!!!!!!
    # also huge reduction in memory, from 16GB to 9.6 !!!!
    # this section is still a HUGE bottleneck, needs to be optimized TO THE CORE
    # nodeStack = zeros(MVector{128,UInt32})

    nodeCur = 1
    nodeStack[nodeCur] = 1
    nodeCur += 1

    intersection::ShapeIntersection = ShapeIntersection(false)

    rayDInv::SVec3f = 1.0f0 ./ ray.direction
    rayDSign::SVec3i = rayDInv .< 0

    # walking stack
    while (nodeCur != 1)
        # if nodeCur > maxNodeCur
        #     println("nodecur")
        # end

        # grab node
        nodeCur -= 1
        @inbounds node = bvh.nodes[nodeStack[nodeCur]]

        # intersect bbox
        if !intersectBbox(ray, rayDInv, node.bbox)
            continue
        end

        # intersect node, switching based on node type
        # for each type, iterate over the primitive list
        if node.internal
            # for internal nodes, attempts to proceed along the
            # split axis from smallest to largest nodes
            if (rayDSign[node.axis] != 0)
                nodeStack[nodeCur] = node.start + 0
                nodeCur += 1
                nodeStack[nodeCur] = node.start + 1
                nodeCur += 1
            else
                nodeStack[nodeCur] = node.start + 1
                nodeCur += 1
                nodeStack[nodeCur] = node.start + 0
                nodeCur += 1
            end
        elseif !isempty(shape.triangles)
            for idx = node.start:node.start+node.num-1
                @inbounds pointAindex, pointBindex, pointCindex =
                    @view shape.triangles[bvh.primitives[idx], :]

                @inbounds pointA = SVec3f(
                    shape.positions[pointAindex, 1],
                    shape.positions[pointAindex, 2],
                    shape.positions[pointAindex, 3],
                )
                @inbounds pointB = SVec3f(
                    shape.positions[pointBindex, 1],
                    shape.positions[pointBindex, 2],
                    shape.positions[pointBindex, 3],
                )
                @inbounds pointC = SVec3f(
                    shape.positions[pointCindex, 1],
                    shape.positions[pointCindex, 2],
                    shape.positions[pointCindex, 3],
                )
                pIntersection =
                    intersectPrimitiveTriangle(ray, pointA, pointB, pointC)
                if !pIntersection.hit
                    continue
                end
                intersection = ShapeIntersection(
                    true,
                    bvh.primitives[idx],
                    pIntersection.u,
                    pIntersection.v,
                    pIntersection.distance,
                )
                ray = Ray(
                    ray.origin,
                    ray.direction,
                    ray.tmin,
                    pIntersection.distance,
                )
            end

        elseif !isempty(shape.quads)
            for idx = node.start:node.start+node.num-1
                @inbounds pointAindex, pointBindex, pointCindex, pointDindex =
                    @view shape.quads[bvh.primitives[idx], :]
                @inbounds pointA = SVec3f(
                    shape.positions[pointAindex, 1],
                    shape.positions[pointAindex, 2],
                    shape.positions[pointAindex, 3],
                )
                @inbounds pointB = SVec3f(
                    shape.positions[pointBindex, 1],
                    shape.positions[pointBindex, 2],
                    shape.positions[pointBindex, 3],
                )
                @inbounds pointC = SVec3f(
                    shape.positions[pointCindex, 1],
                    shape.positions[pointCindex, 2],
                    shape.positions[pointCindex, 3],
                )
                @inbounds pointD = SVec3f(
                    shape.positions[pointDindex, 1],
                    shape.positions[pointDindex, 2],
                    shape.positions[pointDindex, 3],
                )

                pIntersection =
                    intersectPrimitiveQuad(ray, pointA, pointB, pointC, pointD)
                if !pIntersection.hit
                    continue
                end
                intersection = ShapeIntersection(
                    true,
                    bvh.primitives[idx],
                    pIntersection.u,
                    pIntersection.v,
                    pIntersection.distance,
                )
                ray = Ray(
                    ray.origin,
                    ray.direction,
                    ray.tmin,
                    pIntersection.distance,
                )
            end
        end

        if findAny && intersection.hit
            return intersection
        end
    end
    return intersection
end

function intersectPrimitiveTriangle(
    ray::Ray,
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
)::PrimitiveIntersection
    edge1 = p1 - p0
    edge2 = p2 - p0
    pvec = cross(ray.direction, edge2)
    det = dot(edge1, pvec)

    if det == 0
        return PrimitiveIntersection(false)
    end

    inverseDet = 1.0f0 / det

    tvec = ray.origin - p0
    u = dot(tvec, pvec) * inverseDet
    if u < 0 || u > 1
        return PrimitiveIntersection(false)
    end

    qvec = cross(tvec, edge1)
    v = dot(ray.direction, qvec) * inverseDet
    if v < 0 || u + v > 1
        return PrimitiveIntersection(false)
    end

    t = dot(edge2, qvec) * inverseDet
    if t < ray.tmin || t > ray.tmax
        return PrimitiveIntersection(false)
    end

    return PrimitiveIntersection(true, u, v, t)
end

function intersectPrimitiveQuad(
    ray::Ray,
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    p3::SVec3f,
)::PrimitiveIntersection
    if (p2 == p3)
        return intersectPrimitiveTriangle(ray, p0, p1, p3)
    end

    isec1 = intersectPrimitiveTriangle(ray, p0, p1, p3)

    isec2 = intersectPrimitiveTriangle(ray, p2, p3, p1)
    if (isec2.hit)
        isec2 = PrimitiveIntersection(
            true,
            1 - isec2.u,
            1 - isec2.v,
            isec2.distance,
        )
    end

    if (isec1.hit && !isec2.hit)
        return isec1
    elseif (isec2.hit && !isec1.hit)
        return isec2
    elseif (isec1.hit && isec2.hit)
        return isec1.distance < isec2.distance ? isec1 : isec2
    else
        return isec1
    end
end

# Intersect a ray with a axis-aligned bounding box
@inline function intersectBbox(ray::Ray, rayDInv::SVec3f, bbox::Bbox3f)::Bool
    itMin::SVec3f = (bbox.min - ray.origin) .* rayDInv
    itMax::SVec3f = (bbox.max - ray.origin) .* rayDInv
    # PROF: this is a decent bottleneck. any way to make this min and max faster?
    # we see that of tmin we only need the max, and of tmax we only need the min
    # tmin::SVec3f = min.(itMin, itMax)
    # tmax::SVec3f = max.(itMin, itMax)
    # t0::Float32 = max(maximum(tmin), ray.tmin)
    # t1::Float32 = min(minimum(tmax), ray.tmax)
    # maybe better version ? 
    # meh doesn't look like much, maybe a tiny bit: 2.8 --> 2.77 seconds
    maxTmin::Float32 = fastMaximum(fastMin.(itMin, itMax))
    minTmax::Float32 = fastMinimum(fastMax.(itMin, itMax))
    t0::Float32 = fastMax(maxTmin, ray.tmin)
    t1::Float32 = fastMin(minTmax, ray.tmax)
    # mabe better to avoid materialization of the broadcast? nah, compiler probably gets it 
    # maxTmin::Float32 = minMax(itMin, itMax)
    # minTmax::Float32 = maxMin(itMin, itMax)
    # t0::Float32 = max(maxTmin, ray.tmin)
    # t1::Float32 = min(minTmax, ray.tmax)
    t1 *= 1.00000024f0 # for double: 1.0000000000000004
    return t0 <= t1
end

@inbounds function fastMinimum(a::SVec3f)
    ifelse(
        a[1] < a[2],
        ifelse(a[1] < a[3], a[1], a[3]),
        ifelse(a[2] < a[3], a[2], a[3]),
    )
    # if a[1] < a[2]
    #     ifelse(a[1] < a[3], a[1], a[3])
    # else
    #     ifelse(a[2] < a[3], a[2], a[3])
    # end
end
@inbounds function fastMaximum(a::SVec3f)
    ifelse(
        a[1] > a[2],
        ifelse(a[1] > a[3], a[1], a[3]),
        ifelse(a[2] > a[3], a[2], a[3]),
    )
    # if a[1] > a[2]
    #     ifelse(a[1] > a[3], a[1], a[3])
    # else
    #     ifelse(a[2] > a[3], a[2], a[3])
    # end
end
@inline function fastMin(a::Float32, b::Float32)
    a_b = a - b
    (signbit(a_b) || isnan(a)) ? a : b
end
@inline function fastMax(a::Float32, b::Float32)
    b_a = b - a
    (signbit(b_a) || isnan(a)) ? a : b
end

function intersectScene(ray::Ray, scene::Scene)::Intersection

    # in the future this will be a BVH
    intersection = Intersection(false)
    for (instanceIndex, instance) in enumerate(scene.instances)
        shape = scene.shapes[instance.shapeIndex]
        for (triangleIndex, (pointAindex, pointBindex, pointCindex)) in
            enumerate(eachcol(transpose(shape.triangles)))
            @inbounds pointA = SVec3f(
                shape.positions[pointAindex, 1],
                shape.positions[pointAindex, 2],
                shape.positions[pointAindex, 3],
            )
            @inbounds pointB = SVec3f(
                shape.positions[pointBindex, 1],
                shape.positions[pointBindex, 2],
                shape.positions[pointBindex, 3],
            )
            @inbounds pointC = SVec3f(
                shape.positions[pointCindex, 1],
                shape.positions[pointCindex, 2],
                shape.positions[pointCindex, 3],
            )
            triangle = Triangle(
                transformPoint(instance.frame, pointA),
                transformPoint(instance.frame, pointB),
                transformPoint(instance.frame, pointC),
            )

            hit::Intersection =
                intersectTriangle(ray, triangle, instanceIndex, triangleIndex)
            if hit.hit
                ray = Ray(ray.origin, ray.direction, ray.tmin, hit.distance)
                intersection = hit
            end
        end

        for (quadIndex, (pointAindex, pointBindex, pointCindex, pointDindex)) in
            enumerate(eachcol(transpose(shape.quads)))
            @inbounds pointA = SVec3f(
                shape.positions[pointAindex, 1],
                shape.positions[pointAindex, 2],
                shape.positions[pointAindex, 3],
            )
            @inbounds pointB = SVec3f(
                shape.positions[pointBindex, 1],
                shape.positions[pointBindex, 2],
                shape.positions[pointBindex, 3],
            )
            @inbounds pointC = SVec3f(
                shape.positions[pointCindex, 1],
                shape.positions[pointCindex, 2],
                shape.positions[pointCindex, 3],
            )
            @inbounds pointD = SVec3f(
                shape.positions[pointDindex, 1],
                shape.positions[pointDindex, 2],
                shape.positions[pointDindex, 3],
            )
            quad = Quad(
                transformPoint(instance.frame, pointA),
                transformPoint(instance.frame, pointB),
                transformPoint(instance.frame, pointC),
                transformPoint(instance.frame, pointD),
            )

            hit::Intersection =
                intersectQuad(ray, quad, instanceIndex, quadIndex)
            if hit.hit
                ray = Ray(ray.origin, ray.direction, ray.tmin, hit.distance)
                intersection = hit
            end
        end
    end
    return intersection
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
    if textureIdx == -1
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
        s = clamp(textureX, 0, 1) * sizeX
        t = clamp(textureY, 0, 1) * sizeY
    else
        s = rem(textureX, 1) * sizeX
        if (s <= 0)
            s += sizeX
        end
        t = rem(textureY, 1) * sizeY
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

function sampleCamera(
    camera::Camera,
    i::Int64,
    j::Int64,
    imwidth::Int64,
    imheight::Int64,
)
    # f1, f2 = rand(Float32, 2)
    # PROF: two separate randoms are much better, we don't create a vector 
    f1 = rand(Float32)
    f2 = rand(Float32)

    u::Float32 = (i + f1) / imwidth
    v::Float32 = (j + f2) / imheight
    #u = clamp(u, 0, 1)
    #v = clamp(v, 0, 1)

    ray = evalCamera(camera, u, v)
    return ray
end

# takes a camera, image coordinates (uv) and generates a ray connecting them
# PROF: this function takes 2166 ms on profview run(1920, 1080, 2)
# PROF: after fixing takes 437 ms 
function evalCamera(camera::Camera, u::Float32, v::Float32)
    film::SVec2f = (
        camera.aspect >= 1 ?
        SVec2f(camera.film, camera.film / camera.aspect) :
        SVec2f(camera.film * camera.aspect, camera.film)
    )

    # PROF: AYAYAYAYAY once again the square brackets no no no
    # q = SVec3f([film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens])
    q = SVec3f(film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens)

    # ray direction through the lens center
    # FIX: normalize is from LinearAlgebra, we need our Norm
    # i think adding the explicit module is the easiest solution to understand
    dc::SVec3f = -Algebra.norm(q)

    # generate random lens_uv
    # uLens, vLens = rand(Float32, 2)
    # PROF: 2 separate rands are better
    uLens = rand(Float32)
    vLens = rand(Float32)

    uLens, vLens = sampleDisk(uLens, vLens)

    # point on the lens
    # PROF: AYAYAYAYAYAYAAYAYAYAYAAY look at those square brackets no no no no 
    # we are making a heap dynamic vector and trashing it right away no no no
    # pointOnLens = SVec3f([
    #     uLens * camera.aperture / 2.0,
    #     vLens * camera.aperture / 2.0,
    #     0,
    # ])
    # corrected version
    pointOnLens =
        SVec3f(uLens * camera.aperture / 2.0, vLens * camera.aperture / 2.0, 0)

    # # PROF: lets destructure this point on lens
    # pointOnLensX = uLens * camera.aperture / 2.0f0
    # pointOnLensY = vLens * camera.aperture / 2.0f0
    # # pointOnLensZ = 0.0f0

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    # FIX: problema: la nostra funzione si chiama norm. questa è normalize di 
    # LinearAlgebra. Non sono la stessa cosa. Come fa a funzionare?
    # ok, la norm di pellacini (e nostra) fa la l-infinity norm, quella di 
    # linearalgebra fa di default la l-2
    direction = Algebra.norm(pointOnFocusPlane - pointOnLens)
    # PROF: destructure this
    # direction = normalize(
    #     SVec3f(
    #         pointOnFocusPlane[1] - pointOnLensX,
    #         pointOnFocusPlane[2] - pointOnLensY,
    #         pointOnFocusPlane[3],
    #     ),
    # )

    ray_origin = transformPoint(camera.frame, pointOnLens)
    ray_direction = transformDirection(camera.frame, direction)

    return Ray(ray_origin, ray_direction)
end
