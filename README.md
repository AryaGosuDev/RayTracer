# RayTracer

This is a custom coded c++ project that is an Offline Photo Realistic Ray Tracer program with advanced rendering effects. Absolutely no third party code / libs were used in this project.

::: FEATURES :::

1: Custom scene description file reader - Must build a scene using a custom scene description format that will define the universal scene lighting and other important global effects. Each item must also be input as a row input into the file along with its custom attributes inc...color, reflectivity, transparency, softness / roughness for BRDF, refraction, emissivity for physically based rendering.

2: N-order OBJ file reader - able to read objects as an OBJ taking into consideration that the order of faces might be > 3. In fact, they can be faces of > 3 up to any N-order faces.

3: Phong interpolation

4: BSP - implemented a bounding volume hierarchy (BVH) model known as a Binary Space Partition to accelerate ray-face intersection.

5: Ray-Face interaction: Using a Möller-Trumbore algorithm for fast ray-face interception.

6: Advanced Shadow Effects - Shadows, Soft Shadows, AO ( Ambient Occlusion )...

7: Photo-Realistic rendering techniques - Smooth / Rough Surface for BRDF micro faceted surfaces, reflectivity of surfaces, glass / water Snell’s refraction…

8: Global illumination - Implemented a radiosity global illumination model based on a n^2 form factor of the scene geometry along with a Gauss–Seidel matrix problem solver, taking into consideration the physical based characteristics of models by inputting their emissivity/absorptivity values in the scene description file.

9: Tone Mapping

10: Linear algebra techniques : Gauss–Seidel solver, inverse matrix solver, SVD ( singular value decomp )...

Rendering examples of this project can be found in the RayTracingSamples Repository
