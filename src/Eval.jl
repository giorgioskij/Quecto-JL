module Eval
using ..Algebra
using ..World
using ..Intersect
using ..Types
using StaticArrays: dot, cross

export evalNormal,
    evalEnvironment,
    evalTexture,
    evalTexcoord,
    evalShadingPosition,
    evalPosition

function evalNormal(
    shape::Shape,
    intersection::Intersection,
    frame::Frame,
    outgoing::SVec3f,
    materialType::String,
)::SVec3f
    if isempty(shape.normals)
        normal = computeNormal(shape, intersection, frame)
    elseif !isempty(shape.triangles)
        @inbounds t = shape.triangles[intersection.elementIndex]
        @inbounds normal = transformNormal(
            frame,
            norm(
                interpolateTriangle(
                    shape.normals[t.x],
                    shape.normals[t.y],
                    shape.normals[t.z],
                    intersection.u,
                    intersection.v,
                ),
            ),
        )
    elseif !isempty(shape.quads)
        @inbounds q = shape.quads[intersection.elementIndex]
        @inbounds normal = transformNormal(
            frame,
            norm(
                interpolateQuad(
                    shape.normals[q.x],
                    shape.normals[q.y],
                    shape.normals[q.z],
                    shape.normals[q.w],
                    intersection.u,
                    intersection.v,
                ),
            ),
        )
        # TODO normalmap
        # TODO refractive material
        # return ifelse(dot(normal, outgoing) >= 0, normal, -normal)
        # return normal
    else
        error("Only triangles and quads right now")
    end

    if materialType == "refractive"
        return normal
    end

    return dot(normal, outgoing) >= 0 ? normal : -normal
end

function computeNormal(shape::Shape, intersection::Intersection, frame::Frame)
    if !isempty(shape.triangles)
        @inbounds t = shape.triangles[intersection.elementIndex]
        @inbounds return transformNormal(
            frame,
            computeTriangleNormal(
                shape.positions[t.x],
                shape.positions[t.y],
                shape.positions[t.z],
            ),
        )
    elseif !isempty(shape.quads)
        @inbounds q = shape.quads[intersection.elementIndex]
        @inbounds return transformNormal(
            frame,
            computeQuadNormal(
                shape.positions[q.x],
                shape.positions[q.y],
                shape.positions[q.z],
                shape.positions[q.w],
            ),
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
    return norm(cross(pointB - pointA, pointC - pointA))
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

# function evalShadingNormal(
#     scene::Scene,
#     intersection::Intersection,
#     outgoing::SVec3f,
# )::SVec3f
#     instance = scene.instances[intersection.instanceIndex]
#     element = intersection.elementIndex
#     u, v = intersection.u, intersection.v
#     shape = scene.shapes[intersection.shapeIndex]
#     material = scene.materials[instance.materialIndex]
#     if !isempty(shape.triangles) || !isempty(shape.triangles)
#         normal = evalNormal(scene, instance, element, u, v)
#         if material.normalTex != -1
#             normal = evalNormalMap(scene, instance, element, u, v)
#         end
#         if material.type == "refractive"
#             return normal
#         end
#         return ifelse(dot(normal, outgoing) >= 0, normal, -normal)
#         # ignoring lines and points
#     else
#         return SVec3f(0, 0, 0)
#     end
# end

# function evalNormal(
#     scene::Scene,
#     instance::Instance,
#     element::Int,
#     u::Float32,
#     v::Float32,
# )
#     shape = scene.shapes[instance.shapeIndex]
#     if isempty(shape.normals)
#         return computeNormal(shape, element, u, v)
#     end
# end

function evalEnvironment(scene::Scene, direction::SVec3f)::SVec3f
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

    # TEMP: try to invert X and Y
    # textureX = 1 - textureX
    textureY = 1 - textureY

    # FIX: here adding srgb to rgb makes the image more blueish, but not as much as yocto. why?
    textureColor = evalTexture(scene, env.emissionTex, textureX, textureY)
    # textureColor = srgbToRgb(textureColor)
    # @show textureColor
    return env.emission * xyz(textureColor)
end

function evalTexture(
    scene::Scene,
    textureIdx::Int32,
    textureX::Float32,
    textureY::Float32,
    # asLinear::Bool = false,
)::SVec4f
    if textureIdx == -1 || textureIdx == 0
        return SVec4f(1, 1, 1, 1)
    end

    texture = scene.textures[textureIdx]
    return evalTexture(
        texture,
        textureX,
        textureY,
        # asLinear,
        texture.nearest,
        texture.clamp,
    )
end

function evalTexture(
    texture::Texture,
    textureX::Float32,
    textureY::Float32,
    # asLinear::Bool,
    noInterpolation::Bool,
    clampToEdge::Bool,
)::SVec4f
    sizeY, sizeX = size(texture.image)
    # println("sizeY: ", sizeY)
    # println("sizeX: ", sizeX)
    # if !isempty(texture.image)
    #     sizeY, sizeX = size(texture.image)
    # elseif !isempty(texture.hdrImage)
    #     (texture.hdrImage)
    #     sizeY, sizeX = size(texture.hdrImage)
    # else
    #     error("Texture contains no image")
    # end

    # asLinear = false
    # clampToEdge = texture.clamp
    # noInterpolation = texture.nearest
    # noInterpolation = true
    s = 0.0f0
    t = 0.0f0

    if clampToEdge
        s = clamp(textureX, 0, 1.0f0) * sizeX
        t = clamp(textureY, 0, 1.0f0) * sizeY
    else
        # s = rem(textureX, 1.0f0) * sizeX
        s = mod1(textureX, 1.0f0) * sizeX
        # if (s <= 0)
        if (s < 0)
            s += sizeX
        end
        # t = rem(textureY, 1.0f0) * sizeY
        t = mod1(textureY, 1.0f0) * sizeY
        # if (t <= 0)
        if (t < 0)
            t += sizeY
        end
    end

    # i::Int = clamp(Int(floor(s)), 1, sizeX)
    # j::Int = clamp(Int(floor(t)), 1, sizeY)
    i::Int = clamp(floor(Int, s), 0, sizeX - 1)
    j::Int = clamp(floor(Int, t), 0, sizeY - 1)

    # ii::Int = i + 1
    # jj::Int = j + 1
    # # if ii > sizeX
    # if ii >= sizeX
    #     ii -= sizeX
    # end
    # # if jj > sizeY
    # if jj >= sizeY
    #     jj -= sizeY
    # end

    ii = (i + 1) % sizeX
    jj = (j + 1) % sizeY

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

@inline function lookupTexture(texture::Texture, i::Int, j::Int)::SVec4f
    # j indices the ROW of the texture
    # i the column

    sizeY, sizeX = size(texture.image)
    # return texture.image[sizeY-j+1, i]
    # return texture.image[sizeY-j, i+1]
    return texture.image[j+1, i+1]
    # if !isempty(texture.image)
    #     sizeY, sizeX = size(texture.image)
    #     return texture.image[sizeY-j+1, i]
    #     # color = SVec4f(rgba.r, rgba.g, rgba.b, rgba.alpha)
    # elseif !isempty(texture.hdrImage)
    #     sizeY, sizeX = size(texture.hdrImage)
    #     return texture.hdrImage[sizeY-j+1, i]
    #     # color = SVec4f(rgb.r, rgb.g, rgb.b, 1)
    # else
    #     error("Texture contains no image")
    # end

    # if asLinear && !texture.linear
    #     return srgbToRgb(color)
    # end
end

function evalNormalSphere(ray::Ray, sphereCenter::SVec3f)
    # compute coordinates of point hit
    pointHit::SVec3f = ray.origin + ray.tmin * ray.direction
    # compute normal: n = (p-c) / |p-c|
    normal = unitVector(pointHit - sphereCenter)
    return normal
end

function evalTexcoord(
    scene::Scene,
    instance::Instance,
    elementIndex::Int32,
    u::Float32,
    v::Float32,
)
    shape::Shape = scene.shapes[instance.shapeIndex]
    if isempty(shape.textureCoords)
        return u, v
    end
    if !isempty(shape.triangles)
        @inbounds t = shape.triangles[elementIndex]
        @inbounds return interpolateTriangle(
            shape.textureCoords[t.x],
            shape.textureCoords[t.y],
            shape.textureCoords[t.z],
            u,
            v,
        )
    elseif !isempty(shape.quads)
        @inbounds q = shape.quads[elementIndex]
        @inbounds return interpolateQuad(
            shape.textureCoords[q.x],
            shape.textureCoords[q.y],
            shape.textureCoords[q.z],
            shape.textureCoords[q.w],
            u,
            v,
        )
    else
        error("No triangles or quads in this shape")
    end
end

# TODO: remove, its just a wrapper
function evalShadingPosition(scene::Scene, intersection::Intersection)::SVec3f
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
    element::Int32,
    u::Float32,
    v::Float32,
)::SVec3f
    shape = scene.shapes[instance.shapeIndex]
    if !isempty(shape.triangles)
        @inbounds t = shape.triangles[element]
        return transformPoint(
            instance.frame,
            interpolateTriangle(
                shape.positions[t.x],
                shape.positions[t.y],
                shape.positions[t.z],
                u,
                v,
            ),
        )
    elseif !isempty(shape.quads)
        @inbounds q = shape.quads[element]
        return transformPoint(
            instance.frame,
            interpolateQuad(
                shape.positions[q.x],
                shape.positions[q.y],
                shape.positions[q.z],
                shape.positions[q.w],
                u,
                v,
            ),
        )
    else
        error("No triangles or quads in this shape")
    end
end

end