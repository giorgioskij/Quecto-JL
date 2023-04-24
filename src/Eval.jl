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
    evalPosition,
    evalMaterial,
    lookupTexture,
    evalElementNormal
#evalShadingNormal

# function evalShadingNormal(
#     scene::Scene,
#     instance::Instance,
#     elementIndex::Integer,
#     u::Float32,
#     v::Float32,
#     outgoing::SVec3f,
# )
#     shape = scene.shapes[instance.shapeIndex]
#     material = scene.materials[instance.materialIndex]

#     if !isempty(shape.triangles) || !isempty(shape.quads)
#         normal = evalNormal(scene, instance, elementIndex, u, v)
#         if material.normalTex != -1
#             error("no normal textures yet 🤷‍♂️")
#         end
#         if material.type == "refractive"
#             return normal
#         end
#         return dot(normal, outgoing) >= 0 ? normal : -normal
#     elseif !isempty(shape.lines)
#         normal = evalNormal(scene, instance, elementIndex, u, v)
#         return orthonormalize(outgoing, normal)
#     elseif !isempty(shape.points)
#         return outgoing
#     else
#         error("🤷‍♂️")
#     end
# end

# function evalNormal(
#     scene::Scene,
#     instance::Instance,
#     elementIndex::Integer,
#     u::Float32,
#     v::Float32,
# )
#     shape = scene.shapes[instance.shapeIndex]
#     if isempty(shape.normals)
#         return evalElementNormal(scene, instance, elementIndex)
#     end
#     if !isempty(shape.triangles)
#         @inbounds t = shape.triangles[elementIndex]
#         @inbounds normal = transformNormal(
#             instance.frame,
#             norm(
#                 interpolateTriangle(
#                     shape.normals[t.x],
#                     shape.normals[t.y],
#                     shape.normals[t.z],
#                     u,
#                     v,
#                 ),
#             ),
#         )
#         return normal
#     elseif !isempty(shape.quads)
#         @inbounds q = shape.quads[elementIndex]
#         @inbounds normal = transformNormal(
#             instance.frame,
#             norm(
#                 interpolateQuad(
#                     shape.normals[q.x],
#                     shape.normals[q.y],
#                     shape.normals[q.z],
#                     shape.normals[q.w],
#                     u,
#                     v,
#                 ),
#             ),
#         )
#         return normal
#     elseif !isempty(shape.lines)
#         @inbounds l = shape.lines[elementIndex]
#         @inbounds normal = transformNormal(
#             instance.frame,
#             norm(interpolateLine(shape.normals[l.x], shape.normals[l.y], u)),
#         )
#         return normal
#     elseif !isempty(shape.points)
#         @inbounds normal = transformNormal(
#             instance.frame,
#             norm(shape.normals[shape.points[elementIndex]]),
#         )
#         return normal
#     else
#         error("🤷‍♂️")
#     end
# end
# function evalNormal(
#     shape::Shape,
#     intersection::Intersection,
#     frame::Frame,
#     outgoing::SVec3f,
#     materialType::String,
# )::SVec3f

#     # if there are triangles or quads

#     # if there are lines

# end
function evalNormalMap(
    scene::Scene,
    instance::Instance,
    elementIndex::Integer,
    u::Float32,
    v::Float32,
    normal::SVec3f,
)::SVec3f
    shape = scene.shapes[instance.shapeIndex]
    material = scene.materials[instance.materialIndex]

    TextureU, TextureV = evalTexcoord(scene, instance, elementIndex, u, v)
    # this is always true in our case
    # if material.normalTex != -1 && (!isempty(shape.triangles) || !isempty(shape.quads))

    normalTex = scene.textures[material.normalTex]
    normalMap =
        -1.0f0 .+
        2.0f0 * xyz(evalTexture(normalTex, TextureU, TextureV, false, false))
    tu, tv = evalElementTangents(scene, instance, elementIndex)
    x = tu
    y = tv
    z = normal
    x = orthonormalize(x, z)
    y = norm(cross(z, x))
    frame = Frame(x, y, z, zeroSV3f)
    flipV::Bool = dot(frame.y, tv) < 0

    if !flipV
        normalMap = SVec3f(normalMap.x, normalMap.y * -1, normalMap.z)
    end

    normal = transformNormal(frame, normalMap)

    return normal
end

function evalElementTangents(
    scene::Scene,
    instance::Instance,
    elementIndex::Integer,
)::Tuple{SVec3f,SVec3f}
    shape = scene.shapes[instance.shapeIndex]

    if !isempty(shape.triangles) && !isempty(shape.textureCoords)
        t = shape.triangles[elementIndex]
        tu, tv = triangleTangentsFromuv(
            shape.positions[t.x],
            shape.positions[t.y],
            shape.positions[t.z],
            shape.textureCoords[t.x],
            shape.textureCoords[t.y],
            shape.textureCoords[t.z],
        )

        return transformDirection(instance.frame, tu),
        transformDirection(instance.frame, tv)

    elseif !isempty(shape.quads) && !isempty(shape.textureCoords)
        q = shape.quads[elementIndex]
        tu, tv = quadTangentsFromuv(
            shape.positions[q.x],
            shape.positions[q.y],
            shape.positions[q.z],
            shape.positions[q.w],
            shape.textureCoords[q.x],
            shape.textureCoords[q.y],
            shape.textureCoords[q.z],
            shape.textureCoords[q.w],
            SVec2f(0, 0),
        )
        return transformDirection(instance.frame, tu),
        transformDirection(instance.frame, tv)
    else
        error("🤷‍♂️")
    end
end

function triangleTangentsFromuv(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    uv0::SVec2f,
    uv1::SVec2f,
    uv2::SVec2f,
)
    p = p1 - p0
    q = p2 - p0
    s = SVec2f(uv1.x - uv0.x, uv2.x - uv0.x)
    t = SVec2f(uv1.y - uv0.y, uv2.y - uv0.y)
    div = s.x * t.y - s.y * t.x

    if (div != 0)
        tu =
            SVec3f(
                t.y * p.x - t.x * q.x,
                t.y * p.y - t.x * q.y,
                t.y * p.z - t.x * q.z,
            ) / div
        tv =
            SVec3f(
                s.x * q.x - s.y * p.x,
                s.x * q.y - s.y * p.y,
                s.x * q.z - s.y * p.z,
            ) / div
        return tu, tv
    else
        return SVec3f(1, 0, 0), SVec3f(0, 1, 0)
    end
end

function quadTangentsFromuv(
    p0::SVec3f,
    p1::SVec3f,
    p2::SVec3f,
    p3::SVec3f,
    uv0::SVec2f,
    uv1::SVec2f,
    uv2::SVec2f,
    uv3::SVec2f,
    currentUv::SVec2f,
)
    if currentUv.x + currentUv.y <= 1
        return triangleTangentsFromuv(p0, p1, p3, uv0, uv1, uv3)
    else
        return triangleTangentsFromuv(p2, p3, p1, uv2, uv3, uv1)
    end
end

function evalNormal(
    # shape::Shape,
    # intersection::Intersection,
    # frame::Frame,
    # outgoing::SVec3f,
    # materialType::String,
    # normalTexId::Int32,

    scene::Scene,
    instance::Instance,
    elementIndex::Int32,
    u::Float32,
    v::Float32,
    outgoing::SVec3f,
)::SVec3f
    shape = scene.shapes[instance.shapeIndex]
    material = scene.materials[instance.materialIndex]
    frame = instance.frame

    if !isempty(shape.triangles) || !isempty(shape.quads)
        # if material.normalTex <= 0
        # eval normal
        if isempty(shape.normals)
            return computeNormal(shape, elementIndex, frame)
        elseif !isempty(shape.triangles)
            @inbounds t = shape.triangles[elementIndex]
            @inbounds normal = transformNormal(
                frame,
                norm(
                    interpolateTriangle(
                        shape.normals[t.x],
                        shape.normals[t.y],
                        shape.normals[t.z],
                        u,
                        v,
                    ),
                ),
            )
        elseif !isempty(shape.quads)
            @inbounds q = shape.quads[elementIndex]
            @inbounds normal = transformNormal(
                frame,
                norm(
                    interpolateQuad(
                        shape.normals[q.x],
                        shape.normals[q.y],
                        shape.normals[q.z],
                        shape.normals[q.w],
                        u,
                        v,
                    ),
                ),
            )
        end
        # eval normalmap
        if material.normalTex > 0
            normal = evalNormalMap(scene, instance, elementIndex, u, v, normal)
        end
        # check if refractive
        if material.type == "refractive"
            return normal
        end

        # check direction
        return dot(normal, outgoing) >= 0 ? normal : -normal

        #if lines
    elseif !isempty(shape.lines)
        if isempty(shape.normals)
            return computeNormal(shape, elementIndex, frame)
        end
        @inbounds l = shape.lines[elementIndex]
        @inbounds normal = transformNormal(
            frame,
            norm(interpolateLine(shape.normals[l.x], shape.normals[l.y], u)),
        )
        return orthonormalize(outgoing, normal)

        # if points
    elseif !isempty(shape.points)
        return outgoing
    end

    error("🤷‍♂️")
end

# function evalElementNormal(
#     scene::Scene,
#     instance::Instance,
#     elementIndex::Integer,
# )::SVec3f
#     shape = scene.shapes[instance.shapeIndex]
#     frame = instance.frame
#     if !isempty(shape.triangles)
#         @inbounds t = shape.triangles[elementIndex]
#         @inbounds return transformNormal(
#             frame,
#             computeTriangleNormal(
#                 shape.positions[t.x],
#                 shape.positions[t.y],
#                 shape.positions[t.z],
#             ),
#         )
#     elseif !isempty(shape.quads)
#         @inbounds q = shape.quads[elementIndex]
#         @inbounds return transformNormal(
#             frame,
#             computeQuadNormal(
#                 shape.positions[q.x],
#                 shape.positions[q.y],
#                 shape.positions[q.z],
#                 shape.positions[q.w],
#             ),
#         )
#     elseif !isempty(shape.lines)
#         @inbounds l = shape.lines[elementIndex]
#         @inbounds return transformNormal(
#             frame,
#             lineTangent(shape.positions[l.x], shape.positions[l.y]),
#         )
#     elseif !isempty(shape.points)
#         return SVec3f(0, 0, 1)
#     else
#         error("🤷‍♂️")
#     end
# end

function computeNormal(shape::Shape, elementIndex::Int32, frame::Frame)
    if !isempty(shape.triangles)
        @inbounds t = shape.triangles[elementIndex]
        @inbounds return transformNormal(
            frame,
            computeTriangleNormal(
                shape.positions[t.x],
                shape.positions[t.y],
                shape.positions[t.z],
            ),
        )
    elseif !isempty(shape.quads)
        @inbounds q = shape.quads[elementIndex]
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
#         return zeroSV3f
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
    emission = zeroSV3f
    for env in scene.environments
        emission += evalEnvironment(scene, env, direction)
    end
    return emission
end

@fastmath function evalEnvironment(
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
    # textureY = 1 - textureY

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

@fastmath function evalTexture(
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
            muladd.(
                lookupTexture(texture, i, j),
                (1 - u) * (1 - v),
                muladd.(
                    lookupTexture(texture, i, jj),
                    (1 - u) * v,
                    muladd.(
                        lookupTexture(texture, ii, j),
                        u * (1 - v),
                        lookupTexture(texture, ii, jj) * u * v,
                    ),
                ),
            )
        )
        # return (
        #     lookupTexture(texture, i, j) * (1 - u) * (1 - v) +
        #     lookupTexture(texture, i, jj) * (1 - u) * v +
        #     lookupTexture(texture, ii, j) * u * (1 - v) +
        #     lookupTexture(texture, ii, jj) * u * v
        # )
    end
end

@inline function lookupTexture(texture::Texture, i::Int, j::Int)::SVec4f
    # j indices the ROW of the texture
    # i the column

    # sizeY, sizeX = size(texture.image)
    # return texture.image[sizeY-j+1, i]
    # return texture.image[sizeY-j, i+1]
    return @inbounds texture.image[j+1, i+1]
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
    pointHit::SVec3f = muladd.(ray.tmin, ray.direction, ray.origin)
    #pointHit::SVec3f = ray.origin + ray.tmin * ray.direction
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
)::SVec2f
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
    elseif !isempty(shape.lines)
        @inbounds l = shape.lines[elementIndex]
        @inbounds return interpolateLine(
            shape.textureCoords[l.x],
            shape.textureCoords[l.y],
            u,
        )
    elseif !isempty(shape.points)
        return shape.textureCoords[shape.points[elementIndex]]
    else
        error("🤷‍♂️")
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
    elseif !isempty(shape.lines)
        return evalPosition(scene, instance, element, u, v)
    elseif !isemtpy(shape.points)
        return transformPoint(
            instance.frame,
            shape.positions[shape.points[element]],
        )
    else
        error("🤷‍♂️")
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
    elseif !isempty(shape.lines)
        @inbounds l = shape.lines[element]
        return transformPoint(
            instance.frame,
            interpolateLine(shape.positions[l.x], shape.positions[l.y], u),
        )
    else
        error("No triangles or quads or lines in this shape")
    end
end

function evalMaterial(
    scene::Scene,
    instance::Instance,
    intersection::Intersection,
    minRoughness::Float32 = 0.03f0 * 0.03f0,
)::MaterialPoint
    material::Material = scene.materials[instance.materialIndex]

    # eval texture coordinates
    textureX, textureY = evalTexcoord(
        scene,
        instance,
        intersection.elementIndex,
        intersection.u,
        intersection.v,
    )

    # eval material textures
    materialEmissionTex::SVec4f =
        evalTexture(scene, material.emissionTex, textureX, textureY)
    materialColorTex::SVec4f =
        evalTexture(scene, material.colorTex, textureX, textureY)
    # ignore roughness and scattering textures

    # eval material properties
    emission::SVec3f = material.emission * xyz(materialEmissionTex)
    color::SVec3f = material.color * xyz(materialColorTex)
    roughness::Float32 = material.roughness
    # FIX: roughness needs to be squared for some reason
    roughness = roughness * roughness

    # fix roughness
    if material.type == "matte" ||
       #material.type == "gltfpbr" ||
       material.type == "glossy"
        roughness = clamp(roughness, minRoughness, 1.0f0)
    elseif roughness < minRoughness
        roughness = 0.0f0
    end

    opacity::Float32 = material.opacity * materialColorTex[4]

    return MaterialPoint(
        material.type,
        emission,
        color,
        opacity,
        roughness,
        material.ior,
    )
end

function evalElementNormal(
    scene::Scene,
    instance::Instance,
    elementIndex::Int32,
)::SVec3f
    shape = scene.shapes[instance.shapeIndex]
    if !isempty(shape.triangles)
        @inbounds t = shape.triangles[elementIndex]
        return transformNormal(
            instance.frame,
            computeTriangleNormal(
                shape.positions[t.x],
                shape.positions[t.y],
                shape.positions[t.z],
            ),
        )
    else
        return zeroSV3f
    end
end

end