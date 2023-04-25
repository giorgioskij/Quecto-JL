module Lights
using ..World
using ..Types
using ..Algebra
using ..MaterialsFunctions
using ..Eval
using ..Intersect
using ..Bvh

using StaticArrays: dot, cross

export sampleLights, pdfLights, makeTraceLights, Light

# TODO: move this in World
struct Light
    instance::Int32
    environment::Int32
    elementsCdf::Vector{Float32}
end

# TODO: this is slow because we don't know in advance how many lights we have
# try to think if we know it in advance
function makeTraceLights(scene::Scene)::Vector{Light}
    lights = Vector{Light}()
    # instance lights
    for (i, instance) in enumerate(scene.instances)
        material = scene.materials[instance.materialIndex]
        if material.emission == zeroSV3f
            continue
        end
        shape = scene.shapes[instance.shapeIndex]
        if (isempty(shape.triangles) && isempty(shape.quads))
            continue
        end
        elementsCdf = Vector{Float32}(undef, 0)
        if !isempty(shape.triangles)
            elementsCdf = Vector{Float32}(undef, size(shape.triangles, 1))
            for idx in eachindex(elementsCdf)
                t = shape.triangles[idx]
                elementsCdf[idx] = triangleArea(
                    shape.positions[t.x],
                    shape.positions[t.y],
                    shape.positions[t.z],
                )
                if idx != 1
                    elementsCdf[idx] += elementsCdf[idx-1]
                end
            end
        end
        if !isempty(shape.quads)
            elementsCdf = Vector{Float32}(undef, size(shape.quads, 1))
            for idx in eachindex(elementsCdf)
                q = shape.quads[idx]
                elementsCdf[idx] = quadArea(
                    shape.positions[q.x],
                    shape.positions[q.y],
                    shape.positions[q.z],
                    shape.positions[q.w],
                )
                if idx != 1
                    elementsCdf[idx] += elementsCdf[idx-1]
                end
            end
        end
        push!(lights, Light(i, 0, elementsCdf))
    end

    # environment lights
    for (k, environment) in enumerate(scene.environments)
        if environment.emission == zeroSV3f
            continue
        end

        elementsCdf = Vector{Float32}(undef, 0)
        if environment.emissionTex > 0 # != -1
            texture = scene.textures[environment.emissionTex]
            textureWidth, textureHeight = size(texture.image)
            elementsCdf = Vector{Float32}(undef, textureWidth * textureHeight)
            for idx in eachindex(elementsCdf)
                zidx = idx - 1
                i = floor(Int, zidx % textureWidth)
                j = floor(Int, zidx / textureWidth)
                th = (j + 0.5f0) * pi / textureHeight
                value = lookupTexture(texture, i, j)
                elementsCdf[idx] = maximum(value) * sin(th)
                if idx != 1
                    elementsCdf[idx] += elementsCdf[idx-1]
                end
            end
        end
        push!(lights, Light(0, k, elementsCdf))
    end

    return lights
end

function triangleArea(p0::SVec3f, p1::SVec3f, p2::SVec3f)::Float32
    return Algebra.length(cross(p1 - p0, p2 - p0)) / 2.0f0
end

function quadArea(p0::SVec3f, p1::SVec3f, p2::SVec3f, p3::SVec3f)::Float32
    return triangleArea(p0, p1, p3) + triangleArea(p2, p3, p1)
end

@inline function sampleUniform(size::Int)::Int32
    return floor(Int32, (rand(Float32) * size)) + 1
end

@inline function sampleDiscrete(cdf::Vector{Float32})::Int32
    r = clamp(rand(Float32) * cdf[end], 0.0f0, cdf[end] - 0.00001f0)
    idx = searchsortedlast(cdf, r) + 1 # r is surely < cdf[end] so lastindex(cdf) + 1 is never returned 
    # return idx #clamp(idx, 0, size(cdf, 1) - 1)
    return clamp(idx, 1, size(cdf, 1))
end

# Sample a point uniformly on a triangle returning the baricentric coordinates.
@inline function sampleTriangle(
    ru::Float32,
    rv::Float32,
)::Tuple{Float32,Float32}
    sqrtru = sqrt(ru)
    return 1.0f0 - sqrtru, rv * sqrtru
end

@inline function sampleSphere()::SVec3f
    z = 2.0f0 * rand(Float32) - 1.0f0
    r = sqrt(clamp(1.0f0 - z * z, 0.0f0, 1.0f0))
    phi = 2.0f0 * pi * rand(Float32)
    return SVec3f(r * cos(phi), r * sin(phi), z)
end

function sampleLights(
    scene::Scene,
    lights::Vector{Light},
    position::SVec3f,
)::SVec3f
    lightId = sampleUniform(size(lights, 1))
    light = lights[lightId]
    if light.instance > 0
        instance = scene.instances[light.instance]
        shape = scene.shapes[instance.shapeIndex]
        element = sampleDiscrete(light.elementsCdf)
        u, v = rand(Float32), rand(Float32)
        if !isempty(shape.triangles)
            u, v = sampleTriangle(u, v)
        end
        lightPosition = evalPosition(scene, instance, element, u, v)
        return norm(lightPosition - position)
    elseif light.environment > 0
        environment = scene.environments[light.environment]
        if environment.emissionTex > 0
            emissionTex = scene.textures[environment.emissionTex]
            # emissionTexWidth, emissionTexHeight = size(emissionTex.image)
            emissionTexHeight, emissionTexWidth = size(emissionTex.image)
            idx = sampleDiscrete(light.elementsCdf)
            u, v = ((idx % emissionTexWidth) + 0.5f0) / emissionTexWidth,
            ((idx / emissionTexWidth) + 0.5f0) / emissionTexHeight
            # cidx = CartesianIndices((emissionTexHeight, emissionTexWidth))[idx]
            # u = (cidx[2] + 0.5f0) / emissionTexWidth
            # v = (cidx[1] + 0.5f0) / emissionTexHeight
            return transformDirection(
                environment.frame,
                SVec3f(
                    cos(u * 2.0f0 * pi) * sin(v * pi),
                    cos(v * pi),
                    sin(u * 2.0f0 * pi) * sin(v * pi),
                ),
            )
        else
            return sampleSphere() # never done
        end
    else
        return zeroSV3f
    end
end

# TODO: add to Algebra
@inline function distanceSquared(a::SVec3f, b::SVec3f)::Float32
    return dot(a - b, a - b)
end

function pdfLights(
    scene::Scene,
    bvh::SceneBvh,
    lights::Vector{Light},
    position::SVec3f,
    direction::SVec3f,
)
    pdf = 0.0f0
    for light in lights
        if light.instance > 0
            instance = scene.instances[light.instance]
            #check all intersection
            lpdf = 0.0f0
            nextPosition = position
            for bounce = 1:100
                intersection = intersectInstance(
                    bvh,
                    scene,
                    light.instance,
                    Ray(nextPosition, direction),
                    false,
                    convert(UInt8, Threads.threadid()),
                )
                if !intersection.hit
                    break
                end
                lposition = evalPosition(
                    scene,
                    instance,
                    intersection.elementIndex,
                    intersection.u,
                    intersection.v,
                )
                lnormal = evalElementNormal(
                    scene,
                    instance,
                    intersection.elementIndex,
                )
                area = light.elementsCdf[end]
                lpdf +=
                    distanceSquared(lposition, position) /
                    (abs(dot(lnormal, direction)) * area)
                # continue
                nextPosition = lposition + direction * 1.0f-3
            end
            pdf += lpdf
        elseif light.environment > 0
            environment = scene.environments[light.environment]
            if environment.emissionTex > 0
                emissionTex = scene.textures[environment.emissionTex]
                #emissionTexWidth, emissionTexHeight = size(emissionTex.image)
                emissionTexHeight, emissionTexWidth = size(emissionTex.image)
                wl = transformDirection(inverse(environment.frame), direction)
                u = atan(wl[3], wl[1]) / (2.0f0 * pi)
                v = acos(clamp(wl[2], -1.0f0, 1.0f0)) / pi
                if u < 0
                    u += 1.0f0
                end
                i = clamp(
                    floor(Int, u * emissionTexWidth),
                    0,
                    emissionTexWidth - 1,
                )
                j = clamp(
                    floor(Int, v * emissionTexHeight),
                    0,
                    emissionTexHeight - 1,
                )
                prob =
                    sampleDiscretePdf(
                        light.elementsCdf,
                        j * emissionTexWidth + i,
                    ) / light.elementsCdf[end]
                angle =
                    (2.0f0 * pi / emissionTexWidth) *
                    (pi / emissionTexHeight) *
                    sin(pi * (j + 0.5f0) / emissionTexHeight)
                pdf += prob / angle
            else
                pdf += 1.0f0 / (4.0f0 * pi)
            end
        end
    end
    pdf *= sampleUniformPdf(size(lights, 1))
    return pdf
end

@inline function sampleUniformPdf(size::Int)::Float32
    return 1.0f0 / convert(Float32, size)
end

@inline function sampleDiscretePdf(cdf::Vector{Float32}, idx::Int)::Float32
    if idx == 0
        return cdf[1]
    end
    return (cdf[idx+1] - cdf[idx])
end

# end of module Lights
end
