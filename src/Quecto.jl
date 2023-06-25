module Quecto

# GC.enable_logging(true)

include("Types.jl")
include("Algebra.jl")
include("World.jl")
include("Bvh.jl")
include("Intersect.jl")
include("MaterialsFunctions.jl")
include("Pdfs.jl")
include("Bsdf.jl")
include("Eval.jl")
include("Lights.jl")
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
using .Lights

export trace

function shaderCaller(
    shader::String,
    scene::Scene,
    ray::Ray,
    bvh::SceneBvh,
    maxBounces::Int64,
    lights::Vector{Light},
)::SVec3f
    if lowercase(shader) == "volumetric"
        return shadeVolumetric(scene, ray, bvh, maxBounces, lights)
    elseif lowercase(shader) == "normal"
        return shadeNormal(scene, ray, bvh)
    elseif lowercase(shader) == "color"
        return shadeColor(scene, ray, bvh)
    elseif lowercase(shader) == "raytrace"
        return shadeRaytrace(scene, ray, bvh)
    elseif lowercase(shader) == "material"
        return shadeMaterial(scene, ray, bvh, maxBounces, lights)
    elseif lowercase(shader) == "path"
        return shadePath(scene, ray, bvh, maxBounces, lights)
    elseif lowercase(shader) == "eyelight"
        return shadeEyelight(scene, ray, bvh)
    else
        error("No shader named $shader")
    end
end

"""
Main entry point to the program. The parameter shader is a string between:
eyelight, normal, color, raytrace, material, path, volumetric.
"""
function trace(;
    scenePath::String = "03_texture/texture.json",
    shader::String = "volumetric",
    resolution::Integer = 1280,
    samples::Integer = 512,
    filename::String = "jtrace.png",
    multithreaded::Bool = true,
    quiet::Bool = false,
    maxBounces::Integer = 8,
    camera::Integer = 1,
    displaySampleTime::Bool = false,
)
    if !quiet
        println(
            "~~~~~ SHADER $shader, RESOLUTION $resolution, SAMPLES $samples, ",
            "THREADS $(multithreaded ? Threads.nthreads() : 1) ~~~~~",
        )
    end
    # generate scene
    t = @elapsed scene = loadJsonScene(joinpath(baseDir, scenePath))
    displayStat("Loaded $scenePath", t)

    camera = scene.cameras[camera]

    if camera.aspect >= 1
        width = resolution
        height = Int(round(resolution / camera.aspect))
    else
        height = resolution
        width = Int(round(resolution * camera.aspect))
    end

    # build bvh
    t = @elapsed bvh = makeSceneBvh(scene)
    if !quiet
        displayStat("BVH built", t)
    end

    # generate empty starting image
    image = zeros(SVec4f, height, width)

    t = @elapsed lights = makeTraceLights(scene)

    # @show lights
    # return lights

    if !quiet
        displayStat("Lights built", t)
    end

    # shaderFunc::Function = shadeMaterial

    # call the function to trace samples

    # GC.enable(false)
    t = @elapsed traceSamples!(
        shader,
        image,
        scene,
        width,
        height,
        samples,
        bvh,
        camera,
        multithreaded,
        maxBounces,
        lights,
        displaySampleTime,
    )
    # GC.enable(true)

    if !quiet
        displayStat("Image rendered", t)
    end

    outputPath = "out/" * filename
    t = @elapsed saveImage(outputPath, image, multithreaded)
    if !quiet
        displayStat("Image saved at $outputPath", t)
    end

    nothing
end

function saveImage(
    filename::String,
    image::Matrix{SVec4f},
    multithreaded::Bool,
)::Nothing
    pngImage = zeros(RGBA, size(image))
    if multithreaded
        @inbounds Threads.@threads for i = 1:size(image, 1)
            @inbounds Threads.@threads for j = 1:size(image, 2)
                srgb = rgbToSrgb(clamp01nan.(image[i, j]))
                @inbounds pngImage[i, j] =
                    RGBA(srgb[1], srgb[2], srgb[3], srgb[4])
            end
        end
    else
        @inbounds for i = 1:size(image, 1)
            @inbounds for j = 1:size(image, 2)
                srgb = rgbToSrgb(clamp01nan.(image[i, j]))
                @inbounds pngImage[i, j] =
                    RGBA(srgb[1], srgb[2], srgb[3], srgb[4])
            end
        end
    end
    save(filename, pngImage)

    return nothing
end

function traceSamples!(
    shader::String,
    # this is a such pain in the ass the type here is always a function 
    # and each time this is passed to trace sample it is evaluated and
    # after it got a good type different by Any (with Function go to Any idk) 
    image::Matrix{SVec4f},
    scene::Scene,
    imwidth::Int,
    imheight::Int,
    samples::Int,
    bvh::SceneBvh,
    camera::Camera,
    multithreaded::Bool,
    maxBounces::Int64,
    lights::Vector{Light},
    displaySampleTime::Bool = false,
)::Nothing
    for s = 1:samples
        weight::Float32 = 1.0f0 / s
        t = @elapsed begin
            if multithreaded
                Threads.@threads for i = 1:size(image, 2)
                    Threads.@threads for j = 1:size(image, 1)
                        # Threads.@threads for idx in eachindex(image)
                        # i = ceil(Int, idx / imheight)
                        # j = idx % imheight
                        # j = j == 0 ? imheight : j
                        radiance::SVec3f = traceSample(
                            shader,
                            i,
                            j,
                            scene,
                            camera,
                            imwidth,
                            imheight,
                            bvh,
                            maxBounces,
                            lights,
                        )

                        @inbounds image[j, i] = linInterp(
                            image[j, i],
                            SVec4f(radiance.x, radiance.y, radiance.z, 1.0),
                            weight,
                        )
                    end
                end
            else
                for i = 1:size(image)[2]
                    for j = 1:size(image)[1]
                        radiance::SVec3f = traceSample(
                            shader,
                            i,
                            j,
                            scene,
                            camera,
                            imwidth,
                            imheight,
                            bvh,
                            maxBounces,
                            lights,
                        )

                        @inbounds image[j, i] = linInterp(
                            image[j, i],
                            SVec4f(radiance.x, radiance.y, radiance.z, 1.0),
                            weight,
                        )
                    end
                end
            end
        end
        if displaySampleTime
            displayStat("sample $s / $samples", t)
        end
    end

    return nothing
end

@inline function traceSample(
    shader::String,
    i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int,
    bvh::SceneBvh,
    maxBounces::Int64,
    lights::Vector{Light},
)::SVec3f

    # send a ray
    # ray = sampleCamera(camera, i, j, imwidth, imheight) --- wrong
    ray::Ray = sampleCamera(camera, i - 1, j - 1, imwidth, imheight)

    # call the shader
    # radiance::SVec3f = shader(scene, ray, bvh, maxBounces, lights)
    radiance::SVec3f = shaderCaller(shader, scene, ray, bvh, maxBounces, lights)

    return radiance
end

function prettyPrint(c::Dates.CompoundPeriod)::String
    c = canonicalize(c)
    t = Time(0) + c
    hour = Hour(t)
    minute = Minute(t)
    second = Second(t)
    milli = Millisecond(t)
    micro = Microsecond(t)
    timestring = ""
    for (symbol, value) in zip(
        ["hour", "min", "sec", "ms", "μ"],
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

end