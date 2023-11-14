# openmp-raytracing

# Ray Tracing with OpenMP

This is a simple ray tracing program implemented in C, utilizing the OpenMP framework for parallelization. The program generates a ray-traced image and saves it in the PNM (Portable Any Map) format.

## Compilation

To compile the code, use the following command:

```bash
gcc -std=c11 -fopenmp raytrace.c -o raytrace.c.exe
```

## Usage

Run the compiled executable to generate the ray-traced image and save it as "result.pnm" in the current directory.

```bash
./raytrace.c.exe
```

## Code Structure

The code consists of the following components:

- **Intersection and Shading Functions:** The `intersect` function determines intersections with scene objects, and the `shade` function calculates shading for a given ray.

- **Tile Calculation:** The `calc_tile` function computes a portion of the image (a tile) based on the given parameters.

- **Main Function:** The `main` function orchestrates the ray tracing process, including memory allocation, parallelized tile calculation, and image saving.

## Parallelization with OpenMP

The parallelization is achieved using OpenMP directives. Specifically, the `#pragma omp parallel` and `#pragma omp for` directives are used to distribute the workload among multiple threads.

## Performance Measurement

The program includes timing functionality to measure the execution time and calculate the performance in terms of Megapixels per second (MPixels/s).

## Output

The generated image is saved in the PNM format as "result.pnm". You can use various image viewing tools to visualize the result.

## File Descriptions

- **raytrace.c:** The main C source code file.
- **result.pnm:** The output image file in PNM format.
- **Readme.md:** This readme file providing information about the program.

## Dependencies

The code relies on standard C libraries and the OpenMP framework for parallelization.

## Author
Mohammad Shanur Rahman

## License

This code is provided under the [license] (if applicable). Feel free to modify and distribute it as needed, respecting any applicable licenses for external dependencies.
