#import Pkg 
#Pkg.add("JSON")
#Pkg.add("PlyIO")

using Images
using BenchmarkTools
using LinearAlgebra
using .Types

# include("types.jl")
# include("loader.jl")

# struct Sphere
#     center::UInt32 # indexing Scene.points
#     radius::Float32
#     color::SVec3f # temp

#     Sphere(center, radius) = new(center, radius, (0.925, 0.36, 0.38))
# end

# struct SphereScene
#     points::Vector{SVec3f}
#     spheres::Vector{Sphere}
# end

# struct HitObjectSphere
#     hit::Bool
#     sphereHit::Union{Sphere,Nothing}
#     ray::Union{Ray,Nothing}

#     HitObject(hit, sphereHit, ray) = new(hit, sphereHit, ray)
#     HitObject(hit) = new(hit, nothing, nothing)
# end




# main entry point to the program
function run(width, height, numSamples)

    # reads params and initializes stuff

    # generate scene
    scene = loadJsonScene("julia-pathtracer/02_matte/bunny.json")
    println("Scene loaded")

    # generate empty starting image
    image = zeros(SVec3f, height, width)

    # call the function to trace samples
    traceSamples(image, scene, width, height, numSamples)

    # save the resulting image
    rgbImage = zeros(RGB, size(image))
    for i in 1:size(image)[1], j in 1:size(image)[2]
        # rgbImage[i, j] = RGB(image[i, j][1]^0.45, image[i, j][2]^0.45, image[i, j][3]^0.45)
        rgbImage[i, j] = RGB(image[i, j]...)
    end
    save("prova.png", rgbImage)

end

function traceSamples(image, scene, imwidth, imheight, numSamples)
    camera = scene.cameras[1]
    # loop over pixels
    # TODO: add threads
    println("Starting creation of image...")
    for s in 1:numSamples
        println("Sample $s")
        Threads.@threads for i in 1:size(image)[2]
            Threads.@threads for j in 1:size(image)[1] #Threads.@threads
                color = traceSample(i, j, scene, camera, imwidth, imheight)

                weight::Float32 = 1 / s
                image[j, i] = linInterp(image[j, i], color, weight)

                # image[j, i] += color
                # image[i, j, :] += color
            end
        end
    end
end

function traceSample(i::Int,
    j::Int,
    scene::Scene,
    camera::Camera,
    imwidth::Int,
    imheight::Int
)::SVec3f

    # send a ray
    ray = sampleCamera(camera, i, j, imwidth, imheight)

    # call the shader
    radiance = shader(scene, ray)

    return radiance

end

# function intersectEnvironment(ray)::SVec3f
#     # return {1, 1, 1}
#     return SVec3f([1, 1, 1])
# end

function shaderColor(scene::Scene, ray::Ray)::SVec3f
    hit::HitObject = intersectScene(ray, scene)

    if !hit.hit
        radiance = evalEnvironment()
        return radiance
    end

    ray::Ray = hit.ray

    radiance = SVec3f(0.925, 0.36, 0.38)

    return radiance
end

function shaderNormal(scene::Scene, ray::Ray)::SVec3f
    hit::HitObject = intersectScene(ray, scene)

    if !hit.hit
        radiance = evalEnvironment()
        return radiance
    end

    ray::Ray = hit.ray

    # compute normal of the point hit

    # sphere::Sphere = hit.sphereHit
    # sphereCenter::SVec3f = scene.points[sphere.center]
    # normal = evalNormalSphere(ray, sphereCenter)

    radiance = 0.5 .* (normal .+ 1) .* sphere.color

    return radiance
end


function shaderEyelight(scene::Scene, ray::Ray)

    hit::HitObject = intersectScene(ray, scene)

    if !hit.hit
        radiance = evalEnvironment()
        return radiance
    end

    # find normal
    sphereCenter::SVec3f = scene.points[hit.sphereHit.center]
    normal::SVec3f = evalNormalSphere(hit.ray, sphereCenter)

    # opposite of the ray direction
    outgoing = -hit.ray.direction

    # compute radiance
    radiance = abs(dot(normal, outgoing)) * hit.sphereHit.color

    return radiance
end


function intersectTriangle(ray::Ray, triangle::Triangle)

    edge1 = triangle.y - triangle.x
    edge2 = triangle.z - triangle.x
    pvec = cross(ray.direction, edge2)
    det = dot(edge1, pvec)


    if det == 0
        return HitObject(false)
    end

    inverseDet = 1.0f0 / det

    tvec = ray.origin - triangle.x
    u = dot(tvec, pvec) * inverseDet
    if u < 0 || u > 1
        return HitObject(false)
    end

    qvec = cross(tvec, edge1)
    v = dot(ray.direction, qvec) * inverseDet
    if v < 0 || u + v > 1
        return HitObject(false)
    end

    t = dot(edge2, qvec) * inverseDet
    if t < ray.tmin || t > ray.tmax
        return HitObject(false)
    end

    newRay = Ray(ray.origin, ray.direction, ray.tmin, t)

    return HitObject(true, u, v, newRay)

end

function intersectScene(ray::Ray, scene::Scene)::HitObject

    # in the future this will be a BVH
    hitObject = HitObject(false)
    for instance in scene.instances
        shape = scene.shapes[instance.shape+1] # +1 because array in julia starts at 1
        for triangleIndices in eachcol(shape.triangles')
            @inbounds pointA = shape.positions[triangleIndices[1]+1, :]
            @inbounds pointB = shape.positions[triangleIndices[2]+1, :]
            @inbounds pointC = shape.positions[triangleIndices[3]+1, :]

            # nuovo metodo prova (spoiler: inutile, il compilatore lo fa già)
            # @inbounds pointA = shape.positions'[triangleIndices[1]*3+1:triangleIndices[1]*3+3]
            # @inbounds pointB = shape.positions'[triangleIndices[2]*3+1:triangleIndices[2]*3+3]
            # @inbounds pointC = shape.positions'[triangleIndices[3]*3+1:triangleIndices[3]*3+3]
            triangle = Triangle(transformPoint(instance.frame, SVec3f(pointA)),
                transformPoint(instance.frame, SVec3f(pointB)),
                transformPoint(instance.frame, SVec3f(pointC))
            )

            # VETTORI ASTRATTI INFAMI SBOSBOSBSOBSOBSOBFOASD;IFA; SJDGJASDP'
            # GAWK;' DGBLNKH EGFKLJ;ADJKL;FAJKLSDFJK
            # mai più vettori astratti di merda solo CONCRETI oh yes baby
            # T I P I    C O N C R E T I 🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿🤶🏿
            # 🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝🍝
            # è tutto commentato qua sotto perché così ti spiego bene domani
            # il motivo di tutto

            # pointA = SVec3f(
            #     shape.positions[triangleIndices[1]+1, 1],
            #     shape.positions[triangleIndices[1]+1, 2],
            #     shape.positions[triangleIndices[1]+1, 3]
            # )
            # pointB = SVec3f(
            #     shape.positions[triangleIndices[2]+1, 1],
            #     shape.positions[triangleIndices[2]+1, 2],
            #     shape.positions[triangleIndices[2]+1, 3]
            # )
            # pointC = SVec3f(
            #     shape.positions[triangleIndices[3]+1, 1],
            #     shape.positions[triangleIndices[3]+1, 2],
            #     shape.positions[triangleIndices[3]+1, 3]
            # )
            # triangle = Triangle(transformPoint(instance.frame, pointA),
            #     transformPoint(instance.frame, pointB),
            #     transformPoint(instance.frame, pointC)
            # )


            hit = intersectTriangle(ray, triangle)
            if hit.hit
                ray = hit.ray
                hitObject = hit
            end
        end
    end
    return hitObject
end

function evalEnvironment()
    background = SVec3f(0.105, 0.443, 0.90)
    return background
end

function evalNormalSphere(ray::Ray, sphereCenter::SVec3f)
    # compute coordinates of point hit
    pointHit::SVec3f = ray.origin + ray.tmin * ray.direction
    # compute normal: n = (p-c) / |p-c|
    normal = unitVector(pointHit - sphereCenter)
    return normal
end


function sampleCamera(camera::Camera, i::Int, j::Int, imwidth::Int, imheight::Int)

    f1, f2 = rand(Float32, 2)

    u::Float32 = (i + f1) / imwidth
    v::Float32 = (j + f2) / imheight
    #u = clamp(u, 0, 1)
    #v = clamp(v, 0, 1)

    ray = evalCamera(camera, u, v)
    return ray
end

# takes a camera, image coordinates (uv) and generates a ray connecting them
function evalCamera(camera::Camera, u::Float32, v::Float32)
    film::SVec2f = (camera.aspect >= 1
                    ? SVec2f([camera.film, camera.film / camera.aspect])
                    : SVec2f([camera.film * camera.aspect, camera.film]))

    q = SVec3f([film[1] * (0.5f0 - u), film[2] * (v - 0.5f0), camera.lens])

    # ray direction through the lens center
    dc::SVec3f = -normalize(q)

    # generate random lens_uv
    uLens, vLens = rand(Float32, 2)
    uLens, vLens = sampleDisk(uLens, vLens)

    # point on the lens
    pointOnLens = SVec3f([
        uLens * camera.aperture / 2.0,
        vLens * camera.aperture / 2.0, 0
    ])

    # point on the focus plane
    pointOnFocusPlane = dc * camera.focus / abs(dc[3])

    # correct ray direction to account for camera focusing
    direction = normalize(pointOnFocusPlane - pointOnLens)


    ray_origin = transformPoint(camera.frame, pointOnLens)
    ray_direction = transformDirection(camera.frame, direction)

    return Ray(ray_origin, ray_direction)
end

function normalize(v::SVec3f)::SVec3f
    l = length(v)
    return (l != 0) ? v / l : v
end

function length(v::SVec3f)::Float32
    return sqrt(dot(v, v))
end

function transformPoint(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3] + frame.o
end

function transformVector(frame::Frame, v::SVec3f)::SVec3f
    return frame.x * v[1] + frame.y * v[2] + frame.z * v[3]
end

function transformDirection(frame::Frame, v::SVec3f)::SVec3f
    return normalize(transformVector(frame, v))
end

function unitVector(v::SVec3f)::SVec3f
    return v / length(v)
end

function linInterp(a::SVec3f, b::SVec3f, weight::Float32)
    return a * (1 - weight) + b * weight
end

function sampleDisk(u::Float32, v::Float32)::SVec2f
    r = sqrt(v)
    phi = 2 * pi * u
    return cos(phi) * r, sin(phi) * r
end



const shader = shaderColor

# @time run(64, 64, 2)





#TODO: benchmark iteration order
#TODO: try out if image is better srotolata or not
#TODO: check if explicit inlining is possible
#TODO: check argument passing