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
using Dates
using .Types
using .Algebra
using .Bvh
using .World
using .Shading

export trace, t

function prettyPrint(c::Dates.CompoundPeriod)::String
    c = canonicalize(c)
    t = Time(0) + c
    hour = Hour(t)
    minute = Minute(t)
    second = Second(t)
    milli = Millisecond(t)
    micro = Microsecond(t)
    timestring = ""
    @inbounds for (symbol, value) in zip(
        ["hour", "min", "sec", "m", "μ"],
        (hour, minute, second, milli, micro),
    )
        if value.value == 0 && typeof(value) != Second
            continue
        end
        n = lpad(
            value.value,
            typeof(value) == Millisecond || typeof(value) == Microsecond ?
            3 : 2,
            "0",
        )
        timestring *= "$n $symbol, "
    end
    return timestring[1:Base.length(timestring)-2]
end

function displayStat(
    message::String,
    time::Union{Real,Nothing} = nothing,
    maxsize = 40,
)
    if !isnothing(time)
        message *= ":"
    end
    if Base.length(message) > maxsize
        diff = Base.length(message) - maxsize
        message = message[1:end-diff-3] * "..."
    else
        message = message * (" "^(maxsize - Base.length(message)))
    end

    if isnothing(time)
        println(message)
    else
        c = Dates.CompoundPeriod(Microsecond(round(time * 10^6)))
        timestring = prettyPrint(c)
        # timestring = Time(0) + canonicalize(c)
        println("$message $timestring")
    end
end

# main entry point to the program
function trace(;
    scenePath::String = "03_texture/texture.json",
    shader::String = "eyelight",
    width = 1920,
    samples = 2,
    multithreaded::Bool = true,
)
    println(
        "~~~~~ SHADER $shader, WIDTH $width, SAMPLES $samples, ",
        "THREADS $(Threads.nthreads()) ~~~~~",
    )
    # generate scene
    t = @elapsed scene = loadJsonScene(joinpath(baseDir, scenePath))
    displayStat("Loaded $scenePath", t)

    camera = scene.cameras[1]
    height = Int(round(width / camera.aspect))

    # build bvh
    t = @elapsed bvh = makeSceneBvh(scene)
    displayStat("BVH built", t)

    # generate empty starting image
    image = zeros(SVec4f, height, width)
    imageLinear::Bool = true

    if lowercase(shader) == "eyelight"
        shaderFunc = shaderEyelightBsdf
    elseif lowercase(shader) == "normal"
        shaderFunc = shaderNormal
    elseif lowercase(shader) == "naive"
        shaderFunc = shaderIndirectNaive
    else
        error("No shader named $shader")
    end
    # shader = shaderMaterial

    # call the function to trace samples

    t = @elapsed traceSamples(
        shaderFunc,
        image,
        scene,
        width,
        height,
        samples,
        bvh,
        camera,
        multithreaded,
    )
    displayStat("Image rendered", t)

    outputPath = "out/jtrace.png"
    t = @elapsed saveImage(outputPath, image, imageLinear, multithreaded)
    displayStat("Image saved at $outputPath", t)
end

function saveImage(
    filename::String,
    image::Matrix{SVec4f},
    isLinear::Bool,
    multithreaded::Bool,
)
    pngImage = zeros(RGBA, size(image))
    if multithreaded
        @inbounds Threads.@threads for i = 1:size(image, 1)
            @inbounds Threads.@threads for j = 1:size(image, 2)
                if isLinear
                    srgb = rgbToSrgb(clamp01nan.(image[i, j]))
                    pngImage[i, j] = RGBA(srgb[1], srgb[2], srgb[3], srgb[4])
                else
                    rgb = clamp01nan.(image[i, j])
                    pngImage[i, j] = RGBA(rgb[1], rgb[2], rgb[3], rgb[4])
                    error("dont know what to do now")
                end
            end
        end
    else
        @inbounds for i = 1:size(image, 1)
            @inbounds for j = 1:size(image, 2)
                if isLinear
                    srgb = rgbToSrgb(clamp01nan.(image[i, j]))
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
    @inbounds for s = 1:samples
        if multithreaded
            @inbounds Threads.@threads for i = 1:size(image)[2]
                @inbounds Threads.@threads for j = 1:size(image)[1]
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
            @inbounds for i = 1:size(image)[2]
                @inbounds for j = 1:size(image)[1]
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