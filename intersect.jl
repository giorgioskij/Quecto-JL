module Intersect

using StaticArrays
using ..Types
using ..Bvh
using ..Algebra
using ..World
import ..Jtrace

export Intersection, ShapeIntersection, PrimitiveIntersection, intersectScene

const global shapeBvhDepth = 128
const global masterBvhDepth = 128
const global nthreads = Threads.nthreads()
const global masterNodeStack = MVector{masterBvhDepth * nthreads,UInt32}(undef)
const global shapeNodeStack = MVector{shapeBvhDepth * nthreads,UInt32}(undef)

# intersection of a ray with a shape
struct ShapeIntersection
    hit::Bool
    elementIndex::Int64
    u::Float32
    v::Float32
    distance::Float32

    ShapeIntersection(hit::Bool) = new(hit, -1, 0, 0, 0)
    ShapeIntersection(hit, elementIndex, u, v, distance) =
        new(hit, elementIndex, u, v, distance)
end

# intersection of a ray with a primitive element
struct PrimitiveIntersection
    hit::Bool
    u::Float32
    v::Float32
    distance::Float32
    PrimitiveIntersection(hit, u, v, distance) = new(hit, u, v, distance)
    PrimitiveIntersection(hit) = new(hit, 0, 0, typemax(Float32))
end

# intersection of a ray with the scene
struct Intersection
    hit::Bool
    instanceIndex::Int64
    elementIndex::Int64
    u::Float32
    v::Float32
    distance::Float32

    Intersection(hit, instanceIndex, elementIndex, u, v, distance) =
        new(hit, instanceIndex, elementIndex, u, v, distance)
    Intersection(hit) = new(hit, -1, -1, 0, 0, 0)
end

# intersects a scene using a bvh
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

    # nodeStack = MVector{masterBvhDepth,UInt32}(undef)
    # preallocatedNodeStackForMySon = MVector{shapeBvhDepth,UInt32}(undef)

    # PERFORMANCE: for multithreaded global buffer
    threadid::UInt8 = Threads.threadid()
    startindex = (threadid - 1) * masterBvhDepth + 1
    endindex = startindex + masterBvhDepth - 1
    nodeStack = @view masterNodeStack[startindex:endindex]

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

                sIntersection = intersectShapeBvh!(
                    sceneBvh.shapes[instance.shapeIndex],
                    scene.shapes[instance.shapeIndex],
                    invRay,
                    findAny,
                    threadid,
                    # preallocatedNodeStackForMySon,
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

# intersects a shape using a bvh
function intersectShapeBvh!(
    shapeBvh::ShapeBvh,
    shape::Shape,
    ray::Ray,
    findAny::Bool,
    # nodeStack::MVector, # this is passed to avoid multiple initializations
    threadid::UInt8,
)::ShapeIntersection
    bvh = shapeBvh.bvh

    # check empty
    if isempty(bvh.nodes)
        return ShapeIntersection(false)
    end

    # PERFORMANCE: for multithreaded global buffer
    startindex = (threadid - 1) * shapeBvhDepth + 1
    endindex = startindex + shapeBvhDepth - 1
    nodeStack = @view shapeNodeStack[startindex:endindex]

    nodeCur = 1
    nodeStack[nodeCur] = 1
    nodeCur += 1

    intersection::ShapeIntersection = ShapeIntersection(false)

    rayDInv::SVec3f = 1.0f0 ./ ray.direction
    rayDSign::SVec3i = rayDInv .< 0

    # walking stack
    while (nodeCur != 1)

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

# Intersect a ray with a axis-aligned bounding box
@inline function intersectBbox(ray::Ray, rayDInv::SVec3f, bbox::Bbox3f)::Bool
    itMin::SVec3f = (bbox.min - ray.origin) .* rayDInv
    itMax::SVec3f = (bbox.max - ray.origin) .* rayDInv
    maxTmin::Float32 = fastMaximum(fastMin.(itMin, itMax))
    minTmax::Float32 = fastMinimum(fastMax.(itMin, itMax))
    t0::Float32 = fastMax(maxTmin, ray.tmin)
    t1::Float32 = fastMin(minTmax, ray.tmax)
    t1 *= 1.00000024f0 # for double: 1.0000000000000004
    return t0 <= t1
end

@inbounds function fastMinimum(a::SVec3f)
    ifelse(
        a[1] < a[2],
        ifelse(a[1] < a[3], a[1], a[3]),
        ifelse(a[2] < a[3], a[2], a[3]),
    )
end
@inbounds function fastMaximum(a::SVec3f)
    ifelse(
        a[1] > a[2],
        ifelse(a[1] > a[3], a[1], a[3]),
        ifelse(a[2] > a[3], a[2], a[3]),
    )
end
@inline function fastMin(a::Float32, b::Float32)
    a_b = a - b
    (signbit(a_b) || isnan(a)) ? a : b
end
@inline function fastMax(a::Float32, b::Float32)
    b_a = b - a
    (signbit(b_a) || isnan(a)) ? a : b
end

# intersect a primitive triangle, given as a list of 3d points.
function intersectPrimitiveTriangle(
    ray::Ray,
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
)::PrimitiveIntersection
    edge1 = p1 - p0
    edge2 = p2 - p0
    pvec = StaticArrays.cross(ray.direction, edge2)
    det = StaticArrays.dot(edge1, pvec)

    if det == 0
        return PrimitiveIntersection(false)
    end

    inverseDet = 1.0f0 / det

    tvec = ray.origin - p0
    u = StaticArrays.dot(tvec, pvec) * inverseDet
    if u < 0 || u > 1
        return PrimitiveIntersection(false)
    end

    qvec = StaticArrays.cross(tvec, edge1)
    v = StaticArrays.dot(ray.direction, qvec) * inverseDet
    if v < 0 || u + v > 1
        return PrimitiveIntersection(false)
    end

    t = StaticArrays.dot(edge2, qvec) * inverseDet
    if t < ray.tmin || t > ray.tmax
        return PrimitiveIntersection(false)
    end

    return PrimitiveIntersection(true, u, v, t)
end

# intersect a primitive quad, given as a list of 3d points.
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

# ---------------------- END OF BVH RELATED STUFF ----------------------
# --------- the code below is the older naive implementation  ----------

# intersect a scene without a bvh - naive bruteforce method
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

function intersectTriangle(
    ray::Ray,
    triangle::Triangle,
    instanceIndex::Int64,
    triangleIndex::Int64,
)::Intersection
    edge1 = triangle.y - triangle.x
    edge2 = triangle.z - triangle.x
    pvec = StaticArrays.cross(ray.direction, edge2)
    det = StaticArrays.dot(edge1, pvec)

    if det == 0
        return Intersection(false)
    end

    inverseDet = 1.0f0 / det

    tvec = ray.origin - triangle.x
    u = StaticArrays.dot(tvec, pvec) * inverseDet
    if u < 0 || u > 1
        return Intersection(false)
    end

    qvec = StaticArrays.cross(tvec, edge1)
    v = StaticArrays.dot(ray.direction, qvec) * inverseDet
    if v < 0 || u + v > 1
        return Intersection(false)
    end

    t = StaticArrays.dot(edge2, qvec) * inverseDet
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

# end module Intersect 
end