#pragma once

extern float apeture;
extern float focal_length;
extern int depth_of_field_rays;

void addDepthOfField(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, glm::vec3& color);