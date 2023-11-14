#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

// Array representing the scene objects
static double
objs[] = {
    0.,        0.,        -100.5,    10000.,    0.,       0.,        0.,
    0.25,      0.272166,  0.272166,  0.544331,  .027777,  0.643951,  0.172546,
    0.,        .027777,   0.172546,  0.643951,  0.,       .027777,   -0.371785,
    0.099620,  0.544331,  .027777,   -0.471405, 0.471405, 0.,        .027777,
    -0.643951, -0.172546, 0.,        .027777,   0.099620, -0.371785, 0.544331,
    .027777,   -0.172546, -0.643951, 0.,        .027777,  0.471405,  -0.471405,
    0.,        .027777,   4.,        3.,        2.,       1.,        -4.,
    4.,        -3.,       1.,        5.
};

static double
timing()
{
  struct timeval tp;

  if (gettimeofday(&tp, NULL) == -1) {
    fprintf(stderr, "ERROR: gettimeofday failed.\n");
    exit(EXIT_FAILURE);
  }
  
  return tp.tv_sec + tp.tv_usec * 1e-6; 
}

// Function to find intersections with objects in the scene
static double *
intersect(double x, double y, double z, 
          double dx, double dy, double dz,
          double *maxp) 
{
  double *o = objs;
  double *oo = NULL;
  double max = *maxp;

  for (int i = 0; i < 11; i++) {
    double xx = *o++ - x;
    double yy = *o++ - y;
    double zz = *o++ - z;
    double b = xx * dx + yy * dy + zz * dz;
    double t;
    
    if ((t = b * b - xx * xx - yy * yy - zz * zz + *o++) < 0 ||
        (t = b - sqrt(t)) < 1e-6 || t > max)
      continue;
    oo = o - 4;
    max = t;
  }
  *maxp = max;
  return oo;
}

// Function to calculate shading
static double
shade(double x, double y, double z, 
      double dx, double dy, double dz,
      int de) 
{
  double max = 1e6;
  double c = 0.0;

  double *o = intersect(x, y, z, dx, dy, dz, &max);
  if (o == NULL)
    return 0.0;
    
  x += max * dx;
  y += max * dy;
  z += max * dz;
  double nx = x - *o++;
  double ny = y - *o++;
  double nz = z - *o++;
  double r = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= r;
  ny /= r;
  nz /= r;
  double k = nx * dx + ny * dy + nz * dz;
  double rdx = dx - 2 * k * nx;
  double rdy = dy - 2 * k * ny;
  double rdz = dz - 2 * k * nz;
  
  o = objs + 44;
  
  for (int i = 0; i < 3; i++) {
    double ldx = *o++ - x;
    double ldy = *o++ - y;
    double ldz = *o++ - z;
    
    r = sqrt(ldx * ldx + ldy * ldy + ldz * ldz);
    
    ldx /= r;
    ldy /= r;
    ldz /= r;
    
    if (intersect(x, y, z, ldx, ldy, ldz, &r))
      continue;
    if ((r = ldx * nx + ldy * ny + ldz * nz) < 0)
      continue;
    c += r;
    if ((r = rdx * ldx + rdy * ldy + rdz * ldz) > 0)
      c += 2 * pow(r, 15.);
  }
  
  if (de < 10)
    c += .5 * shade(x, y, z, rdx, rdy, rdz, de + 1);
    
  return c;
}
// Calculates a tile in the image
static void
calc_tile(int size, int xstart, int ystart, 
          int tilesize, char *tile) 
{
  int i = 0;

  for (int y = ystart; y < ystart + tilesize; y++) {
    for (int x = xstart; x < xstart + tilesize; x++) {
      double xx = x / (float)(size - 1);
      double yy = 1.0 - y / (float)(size - 1);
      
      double dx = -0.847569 - xx * 1.30741 - yy * 1.19745;
      double dy = -1.98535  + xx * 2.11197 - yy * 0.741279;
      double dz = -2.72303                 + yy * 2.04606;
      
      double r = sqrt(dx * dx + dy * dy + dz * dz);
      int c = (int)(100.0 
              * shade(2.1, 1.3, 1.7, dx / r, dy / r, dz / r, 0));
      
      if (c < 0)    c = 0;
      if (c > 255)  c = 255;
      
      tile[i++] = c;
    }
  }
}

int main(int argc, char **argv)
{
  const int size = 4000;
  const int tilesize = 400;

  (void)argc; (void)argv;
    
  char *picture = (char *)malloc(size * size * sizeof(char));
  if (picture == NULL) {
    fprintf(stderr, "Could not allocate picture memory!\n");
    exit(1);
  }

  for (int i = 0; i < (int)(size * size * sizeof(char)); ++i) {
    picture[i] = 0;
  }

  double time_start = timing();

  /* number of tiles in x and y direction */
  int xtiles = size / tilesize;
  int ytiles = size / tilesize;

  int tile_count = 0;
  
  #pragma omp parallel
  {
    char *tile = (char *)malloc(tilesize * tilesize * sizeof(char));
    if (tile == NULL) {
      fprintf(stderr, "Could not allocate tile memory!\n");
      exit(1);
    }

    #pragma omp for collapse(2) schedule(runtime)
    for (int yc = 0; yc < xtiles; yc++) {
      for (int xc = 0; xc < ytiles; xc++) {
   
        /* calc one tile */
        calc_tile(size, xc * tilesize, yc * tilesize, tilesize, tile);
   
        /* copy to picture buffer */
        for (int i = 0; i < tilesize; i++) {
          int tilebase = yc * tilesize * tilesize * xtiles + xc * tilesize;
          memcpy((void *)(picture + tilebase + i * tilesize * xtiles),
                 (void *)(tile + i * tilesize), tilesize * sizeof(char));
        }
        //tile_count++;
      }
    }
    
    free(tile);
  }
  
  double duration = timing() - time_start;

  fprintf(stderr, "Time: %lf s, Performance: %lf MPixels/s\n", 
          duration,
          (double)size * size / duration / 1e6);

  FILE *fd = fopen("result.pnm", "w");
  if (fd == NULL) {
    fprintf(stderr, "ERROR: could not open result.pnm\n");
    return 1;
  }

  fprintf(fd, "P5\n%d %d\n255\n", size, size);

  for (int count = 0; count < size * size; count++) {
    fputc(picture[count], fd);
  }

  fclose(fd);

  free(picture);

  return 0;
}
