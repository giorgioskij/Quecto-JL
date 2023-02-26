using JSON
using PlyIO
using .Types
using Images

function loadJsonScene(scenePath::String)
    json = JSON.parsefile(scenePath)

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
        defaultTexture = Texture()
        for (i, element) in enumerate(group)
            if !haskey(element, "uri")
                throw(MissingException("uri not present"))
            end
            textureFilenames[i] = element["uri"]
            textures[i] = Texture(
                defaultTexture.image,
                defaultTexture.hdrImage,
                get(element, "linear", defaultTexture.linear),
                get(element, "nearest", defaultTexture.nearest),
                get(element, "clamp", defaultTexture.clamp),
            )
        end
    end

    # Materials
    if haskey(json, "materials")
        group = json["materials"]
        materials = Vector{Material}(undef, size(group, 1))
        # ignore material_names
        defaultMaterial = Material()
        for (i, element) in enumerate(group)
            colorTex = get(element, "color_tex", defaultMaterial.colorTex)
            if colorTex != -1
                colorTex += 1
            end

            emissionTex =
                get(element, "emission_tex", defaultMaterial.emissionTex)
            if emissionTex != -1
                emissionTex += 1
            end

            roughnessTex =
                get(element, "roughness_tex", defaultMaterial.roughnessTex)
            if roughnessTex != -1
                roughnessTex += 1
            end

            scatteringTex =
                get(element, "scatterin_tex", defaultMaterial.scatteringTex)
            if scatteringTex != -1
                scatteringTex += 1
            end

            normalTex = get(element, "normal_tex", defaultMaterial.normalTex)
            if normalTex != -1
                normalTex += 1
            end

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
                emissionTex,
                colorTex,
                roughnessTex,
                scatteringTex,
                normalTex,
            )
            materials[i] = material
        end
    end

    # SHAPES
    # shapeFilenames = Vector{String}(undef, 0)
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

            # +1 because julia arrays starst at 1
            shapeIndex = get(element, "shape", defaultInstance.shapeIndex)
            if shapeIndex != -1
                shapeIndex += 1
            end
            materialIndex =
                get(element, "material", defaultInstance.materialIndex)
            if materialIndex != -1
                materialIndex += 1
            end

            instance = Instance(frame, shapeIndex, materialIndex)
            # ignore params "name", "material"
            instances[i] = instance
        end
    end

    # ENVIRONMENTS
    environments = Vector{Environment}(undef, 0)
    if haskey(json, "environments")
        group = json["environments"]
        environments = Vector{Environment}(undef, size(group, 1))

        defaultEnv = Environment()
        for (i, element) in enumerate(group)
            f = get(element, "frame", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            frame = Frame(f[1:3], f[4:6], f[7:9], f[10:12])

            emissionTex = get(element, "emission_tex", defaultEnv.emissionTex)
            if emissionTex != -1
                emissionTex += 1
            end
            environments[i] = Environment(
                frame,
                get(element, "emission", defaultEnv.emission),
                emissionTex,
            )
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
    for (i, filename) in enumerate(textureFilenames)
        path = (dirname(scenePath), filename) |> joinpath

        # check that extension is png
        extension = path[findlast(==('.'), path)+1:end]
        if lowercase(extension) == "png"
            image::Matrix{RGBA{N0f8}} = loadTexturePng(path)
            hdrImage = Matrix{RGB{N0f16}}(undef, 0, 0)
            textures[i] = Texture(
                image,
                hdrImage,
                false, # png is not linear
                textures[i].nearest,
                textures[i].clamp,
            )
        elseif lowercase(extension) == "hdr"
            hdrImage::Matrix{RGB{N0f16}} = loadTextureHdr(path)
            image = Matrix{RGB{N0f8}}(undef, 0, 0)
            textures[i] = Texture(
                image,
                hdrImage,
                true, # hdr is linear
                textures[i].nearest,
                textures[i].clamp,
            )
        else
            error("only png textures for now!")
        end
    end

    # ignore load subdivs

    scene = Scene(cameras, instances, shapes, textures, materials, environments)

    return scene
end

# loads a texture image in png format
function loadTexturePng(filename::String)::Matrix{RGBA{N0f8}}
    image = load(filename)
    return image
end

# loads a texture image in hdr format
function loadTextureHdr(filename::String)::Matrix{RGB{N0f16}}
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
    if quadCount > 0
        triangles = Matrix{Int64}(undef, 0, 0)
        quads = Matrix{Int64}(undef, quadCount + triCount, 4)
        Threads.@threads for (i, elem) in collect(enumerate(eachrow(faces)))
            elem = elem[1]
            if size(elem)[end] == 4
                quads[i, :] = transpose(elem)
            else # create fake quad from triangle
                quads[i, :] = [elem[1], elem[2], elem[3], elem[3]]
            end
        end
        # +1 because julia arrays start at 1
        quads .+= 1
    elseif triCount > 0
        quads = Matrix{Int64}(undef, 0, 0)
        triangles = Matrix{Int64}(undef, triCount, 3)
        Threads.@threads for (i, elem) in collect(enumerate(eachrow(faces)))
            elem = elem[1]
            triangles[i, :] = transpose(elem)
        end
        # +1 because julia arrays start at 1
        triangles .+= 1
    else
        error("Zero triangles and zero quads. What?")
    end

    # create shape
    return Shape(triangles, quads, positions, normals, textureCoords)
end
