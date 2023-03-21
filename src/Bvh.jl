# oh yeah baby let's make this big bvh boy
# boss level - difficulty hard!
module Bvh
export makeSceneBvh, SceneBvh, ShapeBvh, Bbox3f
using ..Types
using ..World
using ..Algebra

struct Bbox3f
    min::SVec3f
    max::SVec3f

    Bbox3f() = new(
        SVec3f(typemax(Float32), typemax(Float32), typemax(Float32)),
        SVec3f(typemin(Float32), typemin(Float32), typemin(Float32)),
    )
    Bbox3f(min, max) = new(min, max)
end

# PROF: let's try to make this bvh node lighter
# 16 GB of allocations previosly (all Int64)
# It looks a bit faster: 4.11 -> 4.05 seconds, kinda worth
# I wanted to check if it's faster to waste bits to make it aligned (num: Int16)
# Mh, it's pretty comparable. lets use Int16 so that the struct is exactly 256 bits
struct BvhNode
    bbox::Bbox3f   # 6x32 = 192 bits
    start::UInt32  # 32 bits
    num::UInt16    # 16 bits
    axis::UInt8    # 8 bits
    internal::Bool # 8 bits
    # total: 256 bits

    BvhNode() = new(Bbox3f(), 1, 0, 1, false)
    BvhNode(bbox, start, num, axis, internal) =
        new(bbox, start, num, axis, internal)
end

struct BvhTree
    nodes::Vector{BvhNode}
    primitives::Vector{Int}

    BvhTree() = new(Vector{BvhNode}(undef, 0), Vector{Int}(undef, 0))
    BvhTree(nodes, primitives) = new(nodes, primitives)
end

struct ShapeBvh
    bvh::BvhTree
    ShapeBvh() = new(BvhTree())
    ShapeBvh(bvh) = new(bvh)
end

struct SceneBvh
    bvh::BvhTree
    shapes::Vector{ShapeBvh}

    SceneBvh() = new(BvhTree(), Vector{ShapeBvh}(undef, 0))
    SceneBvh(bvh, shapes) = new(bvh, shapes)
end

function makeSceneBvh(scene::Scene)::SceneBvh

    # build bvh for each shape
    bvhForEachShape = Vector{ShapeBvh}(undef, size(scene.shapes, 1))
    @inbounds for (i, shape) in enumerate(scene.shapes)
        bvhForEachShape[i] = makeShapeBvh(shape)
    end

    # instance bounding boxes
    bboxes = Vector{Bbox3f}(undef, size(scene.instances, 1))

    @inbounds for (i, instance) in enumerate(scene.instances)
        bboxes[i] =
            isempty(bvhForEachShape[instance.shapeIndex].bvh.nodes) ? Bbox3f() :
            transformBbox(
                instance.frame,
                bvhForEachShape[instance.shapeIndex].bvh.nodes[1].bbox,
            )
    end

    # build nodes
    masterBvh::BvhTree = makeBvh(bboxes)

    sceneBvh = SceneBvh(masterBvh, bvhForEachShape)
    return sceneBvh
end

function makeShapeBvh(shape::Shape)::ShapeBvh
    # nTriangles = size(shape.triangles, 1)
    # nQuads = size(shape.quads, 1)
    # bboxes = Vector{Bbox3f}(undef, nTriangles + nQuads)
    bboxes = Bbox3f[]

    # bboxes for triangles
    if !isempty(shape.triangles)
        bboxes = Vector{Bbox3f}(undef, size(shape.triangles, 1))

        for (i, t) in enumerate(shape.triangles)
            @inbounds bboxes[i] = triangleBounds(
                shape.positions[t.x],
                shape.positions[t.y],
                shape.positions[t.z],
            )
        end
    elseif !isempty(shape.quads)
        bboxes = Vector{Bbox3f}(undef, size(shape.quads, 1))
        # bboxes for quads

        for (i, q) in enumerate(shape.quads)
            @inbounds bboxes[i] = quadBounds(
                shape.positions[q.x],
                shape.positions[q.y],
                shape.positions[q.z],
                shape.positions[q.w],
            )
        end
    elseif !isempty(shape.lines)
        bboxes = Vector{Bbox3f}(undef, size(shape.lines, 1))
        for (i, l) in enumerate(shape.lines)
            @inbounds bboxes[i] = lineBounds(
                shape.positions[l.x],
                shape.positions[l.y],
                shape.radius[l.x],
                shape.radius[l.y],
            )
        end
    elseif !isemtpy(shape.points)
        bboxes = Vector{Bbox3f}(undef, size(shape.points, 1))
        for (i, p) in enumerate(shape.points)
            @inbounds bboxes[i] =
                pointBounds(shape.positions[p], shape.radius[p])
        end

    else
        error("🤷‍♂️")
    end

    tree::BvhTree = makeBvh(bboxes)

    shapeBvh::ShapeBvh = ShapeBvh(tree)
    return shapeBvh
end

@inline function pointBounds(p::SVec3f, r::Float32)::Bbox3f
    return Bbox3f(min(p .- r, p .+ r), max(p .- r, p .+ r))
end

@inline function lineBounds(
    pointA::SVec3f,
    pointB::SVec3f,
    radiusA::Float32,
    radiusB::Float32,
)::Bbox3f
    return Bbox3f(
        min.(pointA .- radiusA, pointB .- radiusB),
        max(pointA .+ radiusA, pointB .+ radiusB),
    )
end

@inline function triangleBounds(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
)::Bbox3f
    return Bbox3f(min.(pointA, pointB, pointC), max.(pointA, pointB, pointC))
end

@inline function quadBounds(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
    pointD::SVec3f,
)::Bbox3f
    return Bbox3f(
        min.(pointA, pointB, pointC, pointD),
        max.(pointA, pointB, pointC, pointD),
    )
end

function makeBvh(bboxes::Vector{Bbox3f})::BvhTree
    bvhMaxPrims = 4
    nBoxes = size(bboxes, 1)

    # prepare to build nodes
    nodes::Vector{BvhNode} = Vector{BvhNode}(undef, 0)
    sizehint!(nodes, nBoxes * 2)

    # prepare primitives
    primitives::Vector{Int} = Vector{Int}(1:nBoxes)

    # prepare centers
    centers::Vector{SVec3f} = center.(bboxes)

    # push first node onto the stack
    stack::Vector{SVec3i} = Vector{SVec3i}([SVec3i(1, 1, nBoxes)])
    push!(nodes, BvhNode())

    # create nodes until the stack is empty
    while (!isempty(stack))
        # grab node to work on 
        nodeid, startIdx, endIdx = pop!(stack)
        node = nodes[nodeid]

        # compute bounds
        newMin = SVec3f(typemax(Float32), typemax(Float32), typemax(Float32))
        newMax = SVec3f(typemin(Float32), typemin(Float32), typemin(Float32))
        for i = startIdx:endIdx
            newMin, newMax = merge(newMin, newMax, bboxes[primitives[i]])
        end
        newBbox = Bbox3f(newMin, newMax)

        # split into two children
        # FIX: add >= to compensate for different startindex
        if (endIdx - startIdx + 1 > bvhMaxPrims)
            # get split
            mid, axis =
            # splitMiddle(primitives, bboxes, centers, startIdx, endIdx)
                splitMiddle(primitives, centers, startIdx, endIdx)

            # make and internal node
            newInternal = true
            newAxis = axis
            newNum = 2
            newStart = size(nodes, 1) + 1
            push!(nodes, BvhNode())
            push!(nodes, BvhNode())
            push!(stack, SVec3i(newStart + 0, startIdx, mid))
            push!(stack, SVec3i(newStart + 1, mid + 1, endIdx))

            # [start, mid) [mid, end)
            # [start, mid] [mid, end]
        else
            # make leaf node
            newInternal = false
            newNum = endIdx - startIdx + 1
            newStart = startIdx
            newAxis = node.axis
        end
        newNode = BvhNode(newBbox, newStart, newNum, newAxis, newInternal)
        nodes[nodeid] = newNode
    end

    # cleanup
    sizehint!(nodes, size(nodes, 1))

    bvh = BvhTree(nodes, primitives)
    return bvh
end

function splitMiddle(
    primitives::Vector{Int},
    # bboxes::Vector{Bbox3f},
    centers::Vector{SVec3f},
    startIdx::Int32,
    endIdx::Int32,
)
    boxMin = SVec3f(typemax(Float32), typemax(Float32), typemax(Float32))
    boxMax = SVec3f(typemin(Float32), typemin(Float32), typemin(Float32))
    for i = startIdx:endIdx
        boxMin, boxMax = merge(boxMin, boxMax, centers[primitives[i]])
    end
    cbbox = Bbox3f(boxMin, boxMax)
    csize = cbbox.max - cbbox.min

    if csize == zeroSV3f
        return Int(floor((startIdx + endIdx) / 2)), 1
    end

    # split along largest
    axis = 1
    if csize[1] >= csize[2] && csize[1] >= csize[3]
        axis = 1
    end
    if csize[2] >= csize[1] && csize[2] >= csize[3]
        axis = 2
    end
    if csize[3] >= csize[1] && csize[3] >= csize[2]
        axis = 3
    end

    # split the space in the middle along the largest axis
    split = center(cbbox)[axis]
    middle, primitives = partition_centers(
        primitives,
        centers,
        split,
        axis,
        #(primitive) -> centers[primitive][axis] < split, I love this but it is slow!
        startIdx,
        endIdx,
    )

    # if we were not able to split, just break the primitives in half
    if middle == startIdx || middle == endIdx
        middle = div((startIdx + endIdx), 2)
    end

    return middle, axis
end

# partitions the vector based on the value of the function
function partition_centers(
    v::Vector{Int},
    centers::Vector{SVec3f},
    split::Float32,
    axis::Integer,
    lo::Integer,
    hi::Integer,
)
    iswap = lo

    for i = lo:hi
        if (centers[v[i]][axis] < split)
            if iswap != i
                v[iswap], v[i] = v[i], v[iswap]
            end
            iswap += 1
        end
    end

    return iswap, v
end

@inline function merge(minVec::SVec3f, maxVec::SVec3f, b::SVec3f)
    return min.(minVec, b), max.(maxVec, b)
end

@inline function merge(minVec::SVec3f, maxVec::SVec3f, b::Bbox3f)
    return min.(minVec, b.min), max.(maxVec, b.max)
end

@inline center(b::Bbox3f)::SVec3f = (b.min + b.max) / 2

@inline function transformBbox(frame::Frame, bbox::Bbox3f)
    @inbounds corners = [
        SVec3f(bbox.min[1], bbox.min[2], bbox.min[3]),
        SVec3f(bbox.min[1], bbox.min[2], bbox.max[3]),
        SVec3f(bbox.min[1], bbox.max[2], bbox.min[3]),
        SVec3f(bbox.min[1], bbox.max[2], bbox.max[3]),
        SVec3f(bbox.max[1], bbox.min[2], bbox.min[3]),
        SVec3f(bbox.max[1], bbox.min[2], bbox.max[3]),
        SVec3f(bbox.max[1], bbox.max[2], bbox.min[3]),
        SVec3f(bbox.max[1], bbox.max[2], bbox.max[3]),
    ]
    boxMin = SVec3f(typemax(Float32), typemax(Float32), typemax(Float32))
    boxMax = SVec3f(typemin(Float32), typemin(Float32), typemin(Float32))

    for corner in corners
        boxMin, boxMax = merge(boxMin, boxMax, transformPoint(frame, corner))
    end

    return Bbox3f(boxMin, boxMax)
end

end