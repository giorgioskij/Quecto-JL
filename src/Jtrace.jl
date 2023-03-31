module Jtrace

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

export trace

# main entry point to the program
function trace(;
    scenePath::String = "03_texture/texture.json",
    shader::String = "eyelight",
    width = 1920,
    samples = 2,
    multithreaded::Bool = true,
    quiet::Bool = false,
    maxBounces::Int64 = 128,
)::Nothing
    if !quiet
        println(
            "~~~~~ SHADER $shader, WIDTH $width, SAMPLES $samples, ",
            "THREADS $(Threads.nthreads()) ~~~~~",
        )
    end
    # generate scene
    t = @elapsed scene = loadJsonScene(joinpath(baseDir, scenePath))
    displayStat("Loaded $scenePath", t)

    camera = scene.cameras[1]
    height = Int(round(width / camera.aspect))

    # build bvh
    t = @elapsed bvh = makeSceneBvh(scene)
    if !quiet
        displayStat("BVH built", t)
    end

    # generate empty starting image
    image = zeros(SVec4f, height, width)

    shaderFunc::Function = if lowercase(shader) == "eyelight"
        shadeEyelight
    elseif lowercase(shader) == "normal"
        shadeNormal
    elseif lowercase(shader) == "color"
        shadeColor
    elseif lowercase(shader) == "raytrace"
        shadeRaytrace
    elseif lowercase(shader) == "material"
        shadeMaterial
    elseif lowercase(shader) == "path"
        shadePath
    else
        error("No shader named $shader")
    end
    # shaderFunc::Function = shadeMaterial

    # call the function to trace samples

    # GC.enable(false)
    t = @elapsed traceSamples!(
        shaderFunc,
        image,
        scene,
        width,
        height,
        samples,
        bvh,
        camera,
        multithreaded,
        maxBounces,
    )
    # GC.enable(true)

    if !quiet
        displayStat("Image rendered", t)
    end

    outputPath = "out/jtrace.png"
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
    shader, # this is a such pain in the ass the type here is always a function 
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
        #displayStat("sample $s / $samples", t)
    end

    return nothing
end

@inline function traceSample(
    shader::Function,
    i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int,
    bvh::SceneBvh,
    maxBounces::Int64,
)::SVec3f

    # send a ray
    # ray = sampleCamera(camera, i, j, imwidth, imheight) --- wrong
    ray::Ray = sampleCamera(camera, i - 1, j - 1, imwidth, imheight)

    # call the shader
    radiance::SVec3f = shader(scene, ray, bvh, maxBounces)

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