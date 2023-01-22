using StaticArrays, Images
using BenchmarkTools

# const Color = SVector{3,Float32}
# const Point = SVector{3,Float32}
const SVec3f = SVector{3,Float32}

struct Sphere
    center::UInt32
    radius::Float32
end

struct Scene
    points::Vector{SVec3f}
    spheres::Vector{Sphere}
end

# main entry point to the program
function run(height, width)

    # reads params and initializes stuff


    # generate scene
    scene = generateScene()

    # image = zeros(Float32, height, width, 3)
    image = zeros(SVec3f, height, width)

    # call the function to trace samples
    traceSamples(image, scene)

    # save the resulting image

    rgbImage = zeros(RGB, size(image))
    for i in 1:size(image)[1], j in 1:size(image)[2]
        rgbImage[i, j] = RGB(image[i, j][1]^0.45, image[i, j][2]^0.45, image[i, j][3]^0.45)
    end
    # println(image)
    save("prova.png", rgbImage)
    # save("prova.png", image)


end

function generateScene()::Scene
    # generate point vector
    points = Vector{SVec3f}([[64, 64, 1]])

    # generate sphere
    spheres = Vector{Sphere}([Sphere(1, 5.0)])

    return Scene(points, spheres)
end

function traceSamples(image, scene)
    # loop over pixels
    for i in 1:size(image)[1]
        for j in 1:size(image)[2]
            color = traceSample(i, j, scene)
            image[i, j] += color
            # image[i, j, :] += color
        end
    end
end

function traceSample(i, j, scene)
    # intersect scene
    for sphere in scene.spheres

        center = scene.points[sphere.center]
        if (dist2d(SVec3f([i, j, 0]), center) < sphere.radius)
            color = SVec3f([1, 0, 0])
            return color
        else
            color = SVec3f([0, 0, 0])
            return color
        end
    end
    return color
end

function dist2d(p1::SVec3f, p2::SVec3f)
    d = sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
    # println(d)
    return d
end


run(128, 128)




#TODO: benchmark iteration order