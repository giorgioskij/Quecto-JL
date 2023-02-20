using JSON
using PlyIO
using .Types
using Images

function loadJsonScene(filename::String)
    json = JSON.parsefile(filename)

    # if haskey(json, "asset")
    #     element = json["asset"]
    #     copyright = element["copyright"]
    # end

    # CAMERAS
    if haskey(json, "cameras")
        group = json["cameras"]
        cameras = Vector{Camera}(undef, size(group, 1))
        # ignore camera_names
        defaultCamera = Camera()
        for (i, element) in enumerate(group)
            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            camera = Camera(
                frame,
                # get(element, "frame", defaultCamera.frame),
                get(element, "lens", defaultCamera.lens),
                get(element, "film", defaultCamera.film),
                get(element, "aspect", defaultCamera.aspect),
                get(element, "focus", defaultCamera.focus),
                get(element, "aperture", defaultCamera.aperture),
            )
            # ignore params "ortographics", "name", "lookat"
            cameras[i] = camera
        end
    end

    # TEXTURES
    textureFilenames = Vector{String}(undef, 0)
    textures = Vector{Texture}(undef, 0)
    if haskey(json, "textures")
        group = json["textures"]
        textureFilenames = Vector{String}(undef, size(group, 1))
        textures = Vector{Texture}(undef, size(group, 1))
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            textureFilenames[i] = element["uri"]
            linear = haskey(element, "linear") ? element["linear"] : false
            nearest = haskey(element, "nearest") ? element["nearest"] : false
            clamp = haskey(element, "clamp") ? element["clamp"] : false
            textures[i] =
                Texture(Matrix{RGBA{N0f8}}(undef, 0, 0), linear, nearest, clamp)
        end
    end

    # Materials
    if haskey(json, "materials")
        group = json["materials"]
        materials = Vector{Material}(undef, size(group, 1))
        # ignore material_names
        defaultMaterial = Material()
        for (i, element) in enumerate(group)
            material = Material(
                get(element, "type", defaultMaterial.type),
                get(element, "emission", defaultMaterial.emission),
                get(element, "color", defaultMaterial.color),
                get(element, "metallic", defaultMaterial.metallic),
                get(element, "roughness", defaultMaterial.roughness),
                get(element, "ior", defaultMaterial.ior),
                get(element, "trdepth", defaultMaterial.trdepth),
                get(element, "scattering", defaultMaterial.scattering),
                get(element, "scanisotropy", defaultMaterial.scanisotropy),
                get(element, "opacity", defaultMaterial.opacity),
                get(element, "emissionTex", defaultMaterial.emissionTex),
                get(element, "colorTex", defaultMaterial.colorTex),
                get(element, "roughnessTex", defaultMaterial.roughnessTex),
                get(element, "scatteringTex", defaultMaterial.scatteringTex),
                get(element, "normalTex", defaultMaterial.normalTex),
            )
            materials[i] = material
        end
    end

    # SHAPES
    shapeFilenames = Vector{String}(undef, 0)
    if haskey(json, "shapes")
        group = json["shapes"]
        shapeFilenames = Vector{String}(undef, size(group, 1))
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            shapeFilenames[i] = element["uri"]
        end
    end

    # INSTANCES
    if haskey(json, "instances")
        group = json["instances"]
        instances = Vector{Instance}(undef, size(group, 1))

        defaultInstance = Instance()
        for (i, element) in enumerate(group)
            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            instance = Instance(
                frame,
                # +1 because julia arrays starst at 1
                get(element, "shape", defaultInstance.shapeIndex) + 1,
            )
            # ignore params "name", "material"
            instances[i] = instance
        end
    end

    # LOAD RESOURCES

    # load shapes
    shapes = Vector{Shape}(undef, size(shapeFilenames, 1))
    for (i, filename) in enumerate(shapeFilenames)
        path = (dirname(scenePath), filename) |> joinpath
        shape = loadShape(path)
        shapes[i] = shape
    end

    # load textures
    # textures = Vector{Texture}(undef, size(textureFilenames, 1))
    for (i, filename) in enumerate(textureFilenames)
        path = (dirname(scenePath), filename) |> joinpath

        # check that extension is png
        extension = path[findlast(==('.'), path)+1:end]
        if extension != "png"
            error("only png textures for now!")
        end

        image::Matrix{RGBA{N0f8}} = loadTexturePng(path)

        textures[i] = Texture(
            image,
            textures[i].linear,
            textures[i].nearest,
            textures[i].clamp,
        )
    end

    # ignore load subdivs

    scene = Scene(cameras, instances, shapes, textures, materials)

    return scene
end

# loads a texture image. only png for now
function loadTexturePng(filename::String)::Matrix{RGBA{N0f8}}
    image = load(filename)
    return image
end

function loadShape(filename::String)
    ply = load_ply(filename)

    # check properties
    vertexProperties = plyname.(ply["vertex"].properties)

    # load positions
    if "x" in vertexProperties &&
       "y" in vertexProperties &&
       "z" in vertexProperties
        x = ply["vertex"]["x"]
        y = ply["vertex"]["y"]
        z = ply["vertex"]["z"]
        positions = hcat(x, y, z)
    else
        positions = Matrix{Float32}(undef, 0, 0)
    end

    # load normals
    if "nx" in vertexProperties &&
       "ny" in vertexProperties &&
       "nz" in vertexProperties
        nx = ply["vertex"]["nx"]
        ny = ply["vertex"]["ny"]
        nz = ply["vertex"]["nz"]
        normals = hcat(nx, ny, nz)
    else
        normals = Matrix{Float32}(undef, 0, 0)
    end

    # load texture coordinates
    if "u" in vertexProperties && "v" in vertexProperties
        u = ply["vertex"]["u"]
        v = ply["vertex"]["v"]
        textureCoords = hcat(u, v)
    else
        textureCoords = Matrix{Float32}(undef, 0, 0)
    end

    faces = ply["face"]["vertex_indices"]

    # count triangles and quads
    triCount = count(elem -> (size(elem[1])[end] == 3), eachrow(faces))
    quadCount = count(elem -> (size(elem[1])[end] == 4), eachrow(faces))
    if triCount + quadCount != size(faces)[1]
        throw(MissingException("Only implemented triangles and quads"))
    end

    # load triangles and quads
    triangles = Matrix{Int64}(undef, triCount, 3)
    quads = Matrix{Int64}(undef, quadCount, 4)
    triCounter = 1
    quadCounter = 1
    for elem in eachrow(faces)
        elem = elem[1]
        if size(elem)[end] == 4
            quads[quadCounter, :] = transpose(elem)
            quadCounter += 1
        elseif size(elem)[end] == 3
            triangles[triCounter, :] = transpose(elem)
            triCounter += 1
        else
            throw(MissingException("Only implemented triangles and quads"))
        end
    end

    # +1 because julia arrays start at 1
    triangles .+= 1
    quads .+= 1

    # create shape
    return Shape(triangles, quads, positions, normals, textureCoords)
end
