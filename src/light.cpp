#include "light.h"
#include "config.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>

int NUMBER_OF_SHADOW_SAMPLES = 60;

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    // s1------p----s2
    // alpha    1 - alpha
    // use interpolation to determine the color of the point p
	float alpha = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // the point in the segment
    position = segmentLight.endpoint0 + alpha * (segmentLight.endpoint1 - segmentLight.endpoint0);
    color = (1 - alpha) * segmentLight.color0 + alpha * segmentLight.color1;
}


// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{ 
    float alpha = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    float beta = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // the point in the segment
    position = parallelogramLight.v0 + alpha * parallelogramLight.edge01 + beta * parallelogramLight.edge02;
    color = (1.0f - alpha) * (1.0f - beta) * parallelogramLight.color0 + alpha * (1.0f - beta) * parallelogramLight.color1 + (1.0f - alpha) * beta * parallelogramLight.color2 + alpha * beta * parallelogramLight.color3;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise 
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    auto rayPoint = ray.origin + ray.direction * ray.t;
    Ray rayToLight = Ray(rayPoint, samplePos - rayPoint);
    bvh.intersect(rayToLight, hitInfo, features);
    if (rayToLight.t >= 1) {
        rayToLight.t = 1;
        if (DrawShadowRay)
            drawRay(rayToLight, debugColor);
        return 1.0f;
    } else {
        if (DrawShadowRay)
            drawRay(rayToLight, glm::vec3(1, 0, 0));
        return 0.0f;
    }
}

glm::vec3 getLightPosSegmentLight(const SegmentLight light) {
    return light.endpoint0 + 0.5f * (light.endpoint1 - light.endpoint0);
}

glm::vec3 getLightPosParallelogramLight(const ParallelogramLight light) {
    return light.v0 + 0.5f * (light.edge01 + light.edge02);
}

glm::vec3 getLightPosPointLight(const PointLight light) {
    return light.position;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 l = glm::vec3(0);
        // If shading is enabled, compute the contribution from all lights.
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                if (features.enableHardShadow || features.enableSoftShadow) {
                    if (testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo) > 0.5) {
                        l += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                    }
                } else {
                    l += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                }
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 finalColor { 0,0,0 };
                if (features.enableHardShadow || features.enableSoftShadow) {
                    for (int i = 0; i <= NUMBER_OF_SHADOW_SAMPLES; i++) {
                        glm::vec3 posPoint = { 0, 0, 0 };
                        glm::vec3 colPoint = { 0, 0, 0 };

                        // take random point from segment and calculate it's color and position using interpolation
                        sampleSegmentLight(segmentLight, posPoint, colPoint);
                        PointLight pointTransformedInSegLight = PointLight { posPoint, colPoint };

                        if (testVisibilityLightSample(posPoint, colPoint, bvh, features, ray, hitInfo)) {
                            glm::vec3 computedColor = computeShading(pointTransformedInSegLight.position, pointTransformedInSegLight.color, features, ray, hitInfo);
                            finalColor += computedColor;
                        }
                    }
                    l += finalColor / (float)NUMBER_OF_SHADOW_SAMPLES;
                }
                
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 finalColor { 0, 0, 0 };
                if (features.enableHardShadow || features.enableSoftShadow) {

                    for (int i = 0; i <= NUMBER_OF_SHADOW_SAMPLES; i++) {
                        glm::vec3 posPoint = { 0, 0, 0 };
                        glm::vec3 colPoint = { 0, 0, 0 };
                        // take random point from segment and calculate it's color and position using interpolation
                        sampleParallelogramLight(parallelogramLight, posPoint, colPoint);
                        PointLight pointTransformedInSegLight = PointLight { posPoint, colPoint };

                        if (testVisibilityLightSample(pointTransformedInSegLight.position, colPoint, bvh, features, ray, hitInfo)) {
                            glm::vec3 computedColor = computeShading(pointTransformedInSegLight.position, pointTransformedInSegLight.color, features, ray, hitInfo);
                            finalColor += computedColor;
                        }
                    }
                    l += finalColor / (float)NUMBER_OF_SHADOW_SAMPLES;
                }
            }
        }
        l.x = std::min(l.x, 1.0f);
        l.y = std::min(l.y, 1.0f);
        l.z = std::min(l.z, 1.0f);
        return l;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
