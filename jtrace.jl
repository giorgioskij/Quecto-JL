module Jtrace

include("types.jl")
include("algebra.jl")
include("world.jl")
include("bvh.jl")
include("intersect.jl")
include("eval.jl")
include("shading.jl")
include("loader.jl")

const baseDir = dirname(@__FILE__)
cd(baseDir)

using Images
using BenchmarkTools
using StaticArrays
using .Types
using .Algebra
using .Bvh
using .World
using .Shading

export trace

# main entry point to the program
function trace(scenePath::String, width = 1920, numSamples = 2)
    # generate scene
    scene = loadJsonScene(joinpath(baseDir, scenePath))
    # return scene

    camera = scene.cameras[1]
    height = Int(round(width / camera.aspect))

    # build bvh
    bvh = makeSceneBvh(scene)

    # generate empty starting image
    image = zeros(SVec4f, height, width)
    imageLinear::Bool = true

    shader = shaderEyelight
    #shader = shaderMaterial
    #shader = shaderNormal
    # shader = shaderIndirectNaive

    # call the function to trace samples
    traceSamples(shader, image, scene, width, height, numSamples, bvh, camera)

    saveImage("out/jtrace.png", image, imageLinear)

    # save the resulting image
    # rgbImage = zeros(RGB, size(image))
    # Threads.@threads for i = 1:size(image)[1]
    #     Threads.@threads for j = 1:size(image)[2]
    #         rgbImage[i, j] = RGB(image[i, j]...)
    #     end
    # end
    # save("out/prova.png", rgbImage)
end

function saveImage(filename::String, image::Matrix{SVec4f}, isLinear::Bool)
    # for i = 1:size(image, 1)
    # for j = 1:size(image, 2)
    # pngImage = Matrix{RGB{N0f8}}(undef, size(image, 1), size(image, 2))
    pngImage = zeros(RGB, size(image))
    for i = 1:size(image, 1)
        for j = 1:size(image, 2)
            if isLinear
                srgb = clamp01nan.(rgbToSrgb(image[i, j]))
                pngImage[i, j] = RGB(srgb[1], srgb[2], srgb[3])
            else
                error("dont know what to do now")
            end
        end
    end
    save(filename, pngImage)
end

function rgbToSrgb(rgb::SVec4f)::SVec4f
    return SVec4f(rgbToSrgb(rgb.x), rgbToSrgb(rgb.y), rgbToSrgb(rgb.z), rgb[4])
end

function rgbToSrgb(rgb::Float32)::Float32
    return (rgb <= 0.0031308f0) ? (12.92f0 * rgb) :
           ((1 + 0.055f0) * rgb^(1 / 2.4f0) - 0.055f0)
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
                radiance = traceSample(
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
                # clamp.(linInterp(image[j, i], color, weight), 0.0f0, 1.0f0)
                # linInterp(image[j, i], color, weight)
                    linInterp(
                        image[j, i],
                        SVec4f(radiance.x, radiance.y, radiance.z, 1.0),
                        weight,
                    )
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

end