module Eval
using ..Algebra
using ..World
using ..Intersect
using ..Types
using StaticArrays: dot, cross

export evalNormal

function evalNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if isempty(shape.normals)
        return computeNormal(shape, intersection, frame)
    end

    if !isempty(shape.triangles)
        (indexA, indexB, indexC) =
            @view shape.triangles[intersection.elementIndex, :]

        normalA = SVec3f(@view shape.normals[indexA, :])
        normalB = SVec3f(@view shape.normals[indexB, :])
        normalC = SVec3f(@view shape.normals[indexC, :])
        return transformNormal(
            frame,
            norm(
                interpolateTriangle(
                    normalA,
                    normalB,
                    normalC,
                    intersection.u,
                    intersection.v,
                ),
            ),
        )

    elseif !isempty(shape.quads)
        indexA, indexB, indexC, indexD =
            @view shape.quads[intersection.elementIndex, :]

        normalA = SVec3f(@view shape.normals[indexA, :])
        normalB = SVec3f(@view shape.normals[indexB, :])
        normalC = SVec3f(@view shape.normals[indexC, :])
        normalD = SVec3f(@view shape.normals[indexD, :])
        return transformNormal(
            frame,
            norm(
                interpolateQuad(
                    normalA,
                    normalB,
                    normalC,
                    normalD,
                    intersection.u,
                    intersection.v,
                ),
            ),
        )
    else
        error("Only triangles and quads right now")
    end
end

function computeNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if !isempty(shape.triangles)
        (indexA, indexB, indexC) =
            @view shape.triangles[intersection.elementIndex, :]
        pointA = SVec3f(@view shape.positions[indexA, :])
        pointB = SVec3f(@view shape.positions[indexB, :])
        pointC = SVec3f(@view shape.positions[indexC, :])
        return transformNormal(
            frame,
            computeTriangleNormal(pointA, pointB, pointC),
        )
    elseif !isempty(shape.quads)
        indexA, indexB, indexC, indexD =
            @view shape.quads[intersection.elementIndex, :]
        pointA = SVec3f(@view shape.positions[indexA, :])
        pointB = SVec3f(@view shape.positions[indexB, :])
        pointC = SVec3f(@view shape.positions[indexC, :])
        pointD = SVec3f(@view shape.positions[indexD, :])

        return transformNormal(
            frame,
            computeQuadNormal(pointA, pointB, pointC, pointD),
        )
    else
        error("Only triangles and quads right now")
    end
end

# computes the normal of a triangle
@inline function computeTriangleNormal(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
)
    return norm(cross(pointB .- pointA, pointC .- pointA))
end

# computes the normal of a quad 
function computeQuadNormal(
    pointA::SVec3f,
    pointB::SVec3f,
    pointC::SVec3f,
    pointD::SVec3f,
)
    return norm(
        computeTriangleNormal(pointA, pointB, pointD) +
        computeTriangleNormal(pointC, pointD, pointB),
    )
end

end