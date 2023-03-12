#pragma once
#include "common.h"
#include <framework/ray.h>

// Compute the shading at the intersection point using the Phong model.
const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo);

// Given a ray and a normal (in hitInfo), compute the reflected ray in the specular direction (mirror direction).
const Ray computeReflectionRay(Ray ray, HitInfo hitInfo);

glm::vec3 diffuseOnly(const HitInfo hitInfo, const glm::vec3& lPos, const glm::vec3& vPos, const glm::vec3& lightColor, const Features& features, Ray ray);

glm::vec3 phongSpecularOnly(const HitInfo hitInfo, const glm::vec3& lPos, const glm::vec3& vPos, const glm::vec3& cPos, const glm::vec3& lightColor);