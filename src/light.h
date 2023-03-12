#pragma once
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "shading.h"

extern int NUMBER_OF_SHADOW_SAMPLES;
void sampleSegmentLight (const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color);

void sampleParallelogramLight (const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color);

float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

glm::vec3 getLightPosPointLight(PointLight light);

glm::vec3 getLightPosSegmentLight(SegmentLight light);

glm::vec3 getLightPosParallelogramLight(ParallelogramLight light);



