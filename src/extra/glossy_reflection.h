#pragma once

extern int GLOSSY_SAMPLES;
extern float ROUGHNESS;
extern bool DRAW_GLOSSY_CONE;
glm::vec3 getGlossy(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray refletedRay, Ray incomingRay, HitInfo& hitInfo, int depth);