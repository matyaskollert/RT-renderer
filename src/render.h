#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>
#include <vector>
extern int MAX_RAYS_DEPTH;
extern bool SHOW_ONE_RAY;
extern int SHOW_RAY_NUMBER;

// Forward declarations.
struct Scene;
class Screen;
class Trackball;
class BvhInterface;
struct Features;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = 0);

glm::vec3 environment_mapping(glm::vec3 ray);