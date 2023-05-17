#include "raytracer.h"
#include <string.h>
#include <thread>
#include <unistd.h>

#define THREADS 320
std::thread threads[THREADS];

/*
** Checks if the given ray hits a sphere surface and returns.
** Also returns hit data which contains material information.
*/
bool check_sphere_hit(const std::vector<Sphere> &spheres, const Ray &ray,
                      float t_min, float t_max, Hit &hit) {
  Hit closest_hit;
  bool has_hit = false;
  auto closest_hit_distance = t_max;
  Material material;

  for (std::size_t i = 0; i < spheres.size(); i++) {
    const auto &sphere = spheres[i];
    if (sphere_hit(sphere, ray, t_min, closest_hit_distance, closest_hit)) {
      has_hit = true;
      closest_hit_distance = closest_hit.t;
      material = sphere.material;
    }
  }

  if (has_hit) {
    hit = closest_hit;
    hit.material = material;
  }

  return has_hit;
}

/*
** Traces a ray, returns color for the corresponding pixel.
*/
Vector3 trace_ray(const Ray &ray, const std::vector<Sphere> &spheres,
                  int depth) {
  if (depth <= 0) {
    return Vector3(0, 0, 0);
  }

  Hit hit;
  if (check_sphere_hit(spheres, ray, 0.001f, FLT_MAX, hit)) {
    Ray outgoing_ray;
    Vector3 attenuation;

    if (metal_scater(hit.material, ray, hit, attenuation, outgoing_ray)) {
      auto ray_color = trace_ray(outgoing_ray, spheres, depth - 1);
      return Vector3(ray_color.x * attenuation.x, ray_color.y * attenuation.y,
                     ray_color.z * attenuation.z);
    }

    return Vector3(0, 0, 0);
  }

  Vector3 unit_direction = unit_vector(ray.direction);
  auto t = 0.5 * (unit_direction.y + 1.0);
  return Vector3(1.0, 1.0, 1.0) * (1.0 - t) + Vector3(0.5, 0.7, 1.0) * t;
}

void pixel_compute(int samples, int width, int height, int y,
                   const Camera &camera, int depth, Checksum &checksum,
                   const std::vector<Sphere> &spheres) {
  for (int x = 0; x < width; x++) {
    Vector3 pixel_color(0, 0, 0);
    for (int s = 0; s < samples; s++) {
      auto u = (float)(x + random_float()) / (width - 1);
      auto v = (float)(y + random_float()) / (height - 1);
      auto r = get_camera_ray(camera, u, v);
      pixel_color += trace_ray(r, spheres, depth);
    }
    compute_color(checksum, pixel_color, samples);
  }
}

int main(int argc, char **argv) {
  int width = IMAGE_WIDTH;
  int height = IMAGE_HEIGHT;
  int samples = 20;
  int depth = 10;

  // Calculating the aspect ratio and creating the camera for the rendering
  const auto aspect_ratio = (float)width / height;
  Camera camera(Vector3(0, 1, 1), Vector3(0, 0, -1), Vector3(0, 1, 0),
                aspect_ratio, 90, 0.0f, 1.5f);

  std::vector<Sphere> spheres;

  readInput();
  create_random_scene(spheres);

  // checksums for each color individually
  Checksum checksum(0, 0, 0);
  Checksum check_vec[THREADS] = {};

  // Iterate over each pixel and trace a ray to calculate the color.
  // This is done for samples amount of time for each pixel.
  // TODO: Try to parallelize this.
  for (int y = height - 1; y >= THREADS - 1; y -= THREADS) {
    for (int i = 0; i < THREADS; i++) {
      threads[i] = std::thread(pixel_compute, samples, width, height, y - i,
                               std::ref(camera), depth, std::ref(check_vec[i]),
                               std::ref(spheres));
    }
    for (int i = 0; i < THREADS; i++) {
      threads[i].join();
      checksum += check_vec[i];
      check_vec[i] = Checksum(0, 0, 0);
    }
  }
  for (int y = (height % THREADS) - 1; y >= 0; y--) {
    threads[y] =
        std::thread(pixel_compute, samples, width, height, y, std::ref(camera),
                    depth, std::ref(check_vec[y]), std::ref(spheres));
  }
  for (int y = (height % THREADS) - 1; y >= 0; y--) {
    threads[y].join();
    checksum += check_vec[y];
  }

  writeOutput(checksum);

  return 0;
}