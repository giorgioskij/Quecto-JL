# julia-pathtracer
A pathtracer in Julia inspired by Yocto-GL


## Yocto-GL structure
**ytrace.cpp:**  
main function `run`, loops over **samples**, calls `trace_samples`
stores empty image in `state.image`

**yocto_trace.cpp:**
- `trace_samples` loops over pixels, calls `trace_sample`. Writes over image
- `trace_sample` runs the selective shader (recursively). Returns nothing, 
    writes on the image

**ytrace.cpp:**   
now image is written, saves with `save_image`


## Scene Structure
- Vector{Int} points
- Vector{Sphere} spheres
- Vector{Triangles} tri
- Camera


Triangles = Vec3i
Spheres = {Int, Float}

## Schedule
- [x] Initialize flow and save first image 
- [x] Define structures for camera, frame, ray, etc...
- [x] Spheres as {index to point, radius}
- [x] Hit spheres with rays
- [x] Normal shading 
- [x] Antialiasing
- [x] Eyelight shading
- [x] Load triangles
- [ ] Evaluate normals
- [ ] Textures
- [ ] Environment maps
- [ ] Hit other stuff with rays


## Timings 


### Current best results
```
~~~~~ SHADER normal, WIDTH 3000, SAMPLES 100, THREADS 20 ~~~~~
Loaded 03_texture/texture.json:          00 sec, 381 m, 666 μ
BVH built:                               00 sec, 052 m, 833 μ
Image rendered:                          23 sec, 139 m, 236 μ
Image saved at out/jtrace.png:           00 sec, 550 m, 824 μ
```


## Benchmarking
In a starting state, with no computation in the shader, representing the image
with SVectors as colors is slower.


## Long term todo:

Care about RNGs