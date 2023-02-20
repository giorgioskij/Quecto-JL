using Images
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using .Types
using .Algebra
using .Bvh

# main entry point to the program
function run(width, height, numSamples)

    # reads params and initializes stuff

    # generate scene
    scene = loadJsonScene(scenePath)

    # build bvh
    bvh = makeSceneBvh(scene)

    # generate empty starting image
    image = zeros(SVec3f, height, width)

    # call the function to trace samples
    traceSamples(image, scene, width, height, numSamples)

    # save the resulting image
    rgbImage = zeros(RGB, size(image))
    for i = 1:size(image)[1], j = 1:size(image)[2]
        rgbImage[i, j] = RGB(image[i, j]...)
    end
    save("out/prova.png", rgbImage)
end

function traceSamples(image, scene, imwidth, imheight, numSamples)
    camera = scene.cameras[1]
    # loop over pixels
    # TODO: add threads
    # println("Starting creation of image...")
    for s = 1:numSamples
        Threads.@threads for i = 1:size(image)[2]
            Threads.@threads for j = 1:size(image)[1] #Threads.@threads
                color = traceSample(i, j, scene, camera, imwidth, imheight)

                weight::Float32 = 1 / s
                image[j, i] = linInterp(image[j, i], color, weight)
            end
        end
    end
end

function traceSample(
    i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int,
)::SVec3f

    # send a ray
    ray = sampleCamera(camera, i, j, imwidth, imheight)

    # call the shader
    radiance = shader(scene, ray)

    return radiance
end

function shaderColor(scene::Scene, ray::Ray)::SVec3f
    intersection::Intersection = intersectScene(ray, scene)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    instance::Instance = scene.instances[intersection.instanceIndex]
    material::Material = scene.materials[instance.shapeIndex]

    radiance = material.color

    return radiance
end

function shaderNormal(scene::Scene, ray::Ray)::SVec3f
    intersection::Intersection = intersectScene(ray, scene)

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

function shaderEyelight(scene::Scene, ray::Ray)::SVec3f
    intersection::Intersection = intersectScene(ray, scene)

    if !intersection.hit
        radiance = evalEnvironment(scene, ray.direction)
        return radiance
    end

    # compute normal of the point hit
    instance::Instance = scene.instances[intersection.instanceIndex]
    frame::Frame = instance.frame
    shape::Shape = scene.shapes[instance.shapeIndex]
    material::Material = scene.materials[instance.shapeIndex]

    normal = evalNormal(shape, intersection, frame)

    outgoing = -ray.direction

    color = material.color
    # radiance = 0.5 .* (normal .+ 1) .* color
    radiance::SVec3f = abs(dot(normal, outgoing)) .* color

    return radiance
end

function intersectTriangle(
    ray::Ray,
    triangle::Triangle,
    instanceIndex::Int64,
    triangleIndex::Int64,
)
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

    return Intersection(true, instanceIndex, triangleIndex, u, v, t, true)
end

function intersectQuad(
    ray::Ray,
    quad::Quad,
    instanceIndex::Int64,
    quadIndex::Int64,
)
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
            false,
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
            false,
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
# function IntersectQuads(ray::Ray, quad::Quad, instanceIndex::Int64, quadIndex::Int64)::Intersection

#     e10 = quad.b - quad.a
#     e11 = quad.c - quad.b
#     e00 = quad.d - quad.a
#     qn = cross(e10, quad.b - quad.c)
#     q00 = quad.a - ray.origin
#     q10 = quad.b - ray.origin
#     a = dot(cross(q00, ray.direction), e00)
#     c = dot(qn, ray.direction)
#     b = dot(cross(q10, ray.direction), e11)

#     b-= a + c
#     det = b^2 - 4f0 * a * c
#     if det < 0
#         return Intersection(false)
#     end
#     det = sqrt(det)
#     u1 = 0f0
#     u2 = 0f0
#     t = ray.tmax
#     if (c == 0)
#         u1 = -a/b
#         u2 = -1f0
#     else
#         u1 = (-b - copysignf(det, b)) / 2f0
#         u2 = a / u1
#         u1 /= c
#     end

#     if (0 <= u1 && u1 <= 1)
#         pa = linInterp(quad.a, quad.b, u1)
#         pb = linInterp(e00, e11, u1)
#         n = cross(ray.direction, pb)
#         det = dot(n, n)
#         n = cross(n, pa)
#         t1 = dot(n, pb)
#         v1 = dot(n, ray.direction)
#         if (t1 > 0 && 0 <= v1 && v1 <= det)
#             t = t1 / det
#             u = u1
#             v = v1 / det
#         end
#     end

#     if (0 <= u2 && u2 <= 1)
#         pa = linInterp(quad.a, quad.b, u2)
#         pb = linInterp(e00, e11, u2)
#         n = cross(ray.direction, pb)
#         det = dot(n, n)
#         n = cross(n, pa)
#         t2 = dot(n, pb) / det
#         v2 = dot(n, ray.direction)
#         if (0 <= v2 && v2 <= det && t > t2 && t2 > 0)
#             t = t2
#             u = u2
#             v = v2 / det
#         end
#     end

# end

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
    nodeStack = zeros(Int, 128)
    nodeStack[nodeCur] = 1
    nodeCur += 1

    # init intersection
    intersection::Intersection = Intersection(false)

    rayDInv = SVec3f(
        1.0f0 / ray.direction.x,
        1.0f0 / ray.direction.y,
        1.0f0 / ray.direction.z,
    )
    rayDSign = SVec3i(ray_dinv.x < 0, ray_dinv.y < 0, ray_dinv.z < 0)

    # walking stack
    while (nodeCur != 1)

        # grab node
        nodeCur -= 1
        node = masterBvh.nodes[node_stack[nodeCur]]

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
            for idx = node.start:node.start+node.num
                instance = scene.instances[bvh.primitives[idx]]
                invRay = transformRay(inverse(instance.frame, true), ray)

                # to understand this part because we create a lot of objects and then 
                # we possibly rewrite them, we need to understand the performance impact
                # and maybe create the objects only at the end of the loop?
                sIntersection = intersectShapeBvh(
                    sceneBvh.shapes[instance.shape],
                    scene.shapes[instance.shape],
                    invRay,
                    findAny,
                )
                if (!sIntersection.hit)
                    continue
                end
                intersection = Intersection(
                    true,
                    masterBvh.primitives[idx],
                    sIntersection.element,
                    sIntersection.u,
                    sIntersection.v,
                    sIntersectisn.distance,
                    sIntersection.isTriangle,
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

function intersectShapeBvh(
    shapeBvh::ShapeBvh,
    shape::Shape,
    invRay::Ray,
    findAny::Bool,
)::ShapeIntersection
    bvh = shapeBvh.bvh

    # check empty
    if isempty(bvh.nodes)
        return ShapeIntersection(false)
    end

    # node stack
    nodeStack = zeros(Int, 128)
    nodeCur = 1
    nodeStack[nodeCur] = 1
    nodeCur += 1

    intersection::ShapeIntersection = ShapeIntersection(false)

    rayDInv = SVec3f(
        1.0f0 / ray.direction.x,
        1.0f0 / ray.direction.y,
        1.0f0 / ray.direction.z,
    )
    rayDSign = SVec3i(ray_dinv.x < 0, ray_dinv.y < 0, ray_dinv.z < 0)

    # walking stack
    while (nodeCur != 1)

        # grab node
        nodeCur -= 1
        node = bvh.nodes[node_stack[nodeCur]]

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
            error("you've gone too far, kid. go back and fix your primitives.")

end

struct ShapeIntersection
    hit::Bool
    elementIndex::Int64
    u::Float32
    v::Float32
    distance::Float32

    isTriangle::Bool

    ShapeIntersection(hit::Bool) = new(hit, -1, 0, 0, 0, true)
    ShapeIntersection(hit, elementIndex, u, v, distance, isTriangle) =
        new(hit, elementIndex, u, v, distance, isTriangle)
end

# Intersect a ray with a axis-aligned bounding box
@inline function intersectBbox(ray::Ray, rayDInv::SVec3f, bbox::Bbox3f)::Bool
    it_min::SVec3f = (bbox.min - ray.o) * ray_dinv
    it_max::SVec3f = (bbox.max - ray.o) * ray_dinv
    tmin::SVec3f = min.(it_min, it_max)
    tmax::SVec3f = max.(it_min, it_max)
    t0::Float32 = max(maximum(tmin), ray.tmin)
    t1::Float32 = min(minimum(tmax), ray.tmax)
    t1 *= 1.00000024f0 # for double: 1.0000000000000004
    return t0 <= t1
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
    wl::SVec3f = transformDirection(inverse(environment.frame), direction)
    textureX = atan2(wl[3], wl[1]) / (2 * pi)
    textureY = acos(clamp(wl[2], -1, 1)) / pi

    if textureX < 0
        textureX += 1
    end

    return env.emission *
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
    if isempty(texture.image)
        return SVec4f(0, 0, 0, 0)
    end
    sizeX, sizeY = size(texture.image)

    asLinear = false
    clampToEdge = texture.clamp
    noInterpolation = texture.nearest
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
    jj::Int = (j + 1) % sizey
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
    rgba = texture.image[i, j]
    color = SVec4f(rgba.r, rgba.g, rgba.b, rgba.alpha)
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
    f1, f2 = rand(Float32, 2)

    u::Float32 = (i + f1) / imwidth
    v::Float32 = (j + f2) / imheight
    #u = clamp(u, 0, 1)
    #v = clamp(v, 0, 1)

    ray = evalCamera(camera, u, v)
    return ray
end

# takes a camera, image coordinates (uv) and generates a ray connecting them
function evalCamera(camera::Camera, u::Float32, v::Float32)
    film::SVec2f = (
        camera.aspect >= 1 ?
        SVec2f([camera.film, camera.film / camera.aspect]) :
        SVec2f([camera.film * camera.aspect, camera.film])
    )

    q = SVec3f([film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens])

    # ray direction through the lens center
    dc::SVec3f = -normalize(q)

    # generate random lens_uv
    uLens, vLens = rand(Float32, 2)
    uLens, vLens = sampleDisk(uLens, vLens)

    # point on the lens
    pointOnLens = SVec3f([
        uLens * camera.aperture / 2.0,
        vLens * camera.aperture / 2.0,
        0,
    ])

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    direction = normalize(pointOnFocusPlane - pointOnLens)

    ray_origin = transformPoint(camera.frame, pointOnLens)
    ray_direction = transformDirection(camera.frame, direction)

    return Ray(ray_origin, ray_direction)
end
