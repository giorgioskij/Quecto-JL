module Jtrace

include("Types.jl")
include("Algebra.jl")
include("World.jl")
include("Bvh.jl")
include("Intersect.jl")
include("Eval.jl")
include("Shading.jl")
include("Loader.jl")

const baseDir = joinpath(dirname(@__FILE__), "../")
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
function trace(;
    scenePath::String = "03_texture/texture.json",
    shaderName::String = "eyelight",
    width = 1920,
    samples = 2,
    multithreaded::Bool = true,
)
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

    if lowercase(shaderName) == "eyelight"
        shader = shaderEyelightBsdf
    elseif lowercase(shaderName) == "normal"
        shader = shaderNormal
    elseif lowercase(shaderName) == "naive"
        shader = shaderIndirectNaive
    else
        error("No shader named $shaderName")
    end
    # shader = shaderMaterial

    # call the function to trace samples
    traceSamples(
        shader,
        image,
        scene,
        width,
        height,
        samples,
        bvh,
        camera,
        multithreaded,
    )

    saveImage("out/jtrace.png", image, imageLinear, multithreaded)

    # save the resulting image
    # rgbImage = zeros(RGB, size(image))
    # Threads.@threads for i = 1:size(image)[1]
    #     Threads.@threads for j = 1:size(image)[2]
    #         rgbImage[i, j] = RGB(image[i, j]...)
    #     end
    # end
    # save("out/prova.png", rgbImage)
end

function saveImage(
    filename::String,
    image::Matrix{SVec4f},
    isLinear::Bool,
    multithreaded::Bool,
)
    # for i = 1:size(image, 1)
    # for j = 1:size(image, 2)
    # pngImage = Matrix{RGB{N0f8}}(undef, size(image, 1), size(image, 2))
    pngImage = zeros(RGBA, size(image))
    if multithreaded
        Threads.@threads for i = 1:size(image, 1)
            Threads.@threads for j = 1:size(image, 2)
                if isLinear
                    srgb = clamp01nan.(rgbToSrgb(image[i, j]))
                    pngImage[i, j] = RGBA(srgb[1], srgb[2], srgb[3], srgb[4])
                else
                    rgb = clamp01nan.(image[i, j])
                    pngImage[i, j] = RGBA(rgb[1], rgb[2], rgb[3], rgb[4])
                    error("dont know what to do now")
                end
            end
        end
    else
        for i = 1:size(image, 1)
            for j = 1:size(image, 2)
                if isLinear
                    srgb = clamp01nan.(rgbToSrgb(image[i, j]))
                    pngImage[i, j] = RGBA(srgb[1], srgb[2], srgb[3], srgb[4])
                else
                    rgb = clamp01nan.(image[i, j])
                    pngImage[i, j] = RGBA(rgb[1], rgb[2], rgb[3], rgb[4])
                    error("dont know what to do now")
                end
            end
        end
    end
    save(filename, pngImage)
end

function traceSamples(
    shader,
    image,
    scene,
    imwidth,
    imheight,
    samples,
    bvh,
    camera,
    multithreaded,
)
    for s = 1:samples
        if multithreaded
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
        else
            for i = 1:size(image)[2]
                for j = 1:size(image)[1]
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