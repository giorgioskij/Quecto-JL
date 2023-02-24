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
    image = zeros(SVec3f, height, width)

    #shader = shaderEyelight
    #shader = shaderMaterial
    shader = shaderNormal

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

end