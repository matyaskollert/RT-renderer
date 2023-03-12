#include <scene.h>
#include <render.h>
#include <bvh_interface.h>
#include <glm/gtc/random.hpp>
#include "glossy_reflection.h"
#include <iostream>
#include <draw.h>

int GLOSSY_SAMPLES = 2;
float ROUGHNESS = 1.0f;
bool DRAW_GLOSSY_CONE = true;

glm::vec3 getGlossy(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray refletedRay, Ray incomingRay, HitInfo& hitInfo, int depth) {
    float a = 1.0f / hitInfo.material.shininess * ROUGHNESS;
    glm::vec3 U = glm::normalize(glm::cross(refletedRay.direction, incomingRay.direction));
    glm::vec3 V = glm::normalize(glm::cross(refletedRay.direction, U));
    glm::vec3 finalColor(0.0f);

    int realised = 0;
    for (int i = 0; i < GLOSSY_SAMPLES; i++, realised++) {
        float u = ( - a) / 2 + glm::linearRand(0.0f, 1.0f) * a;
        float v = ( - a) / 2 + glm::linearRand(0.0f, 1.0f) * a;
        
        auto newNormal = glm::normalize(refletedRay.direction + u * U + v * V);

        // Checks if the angle would make a ray go behind an object.
        if (glm::dot(newNormal, hitInfo.normal) <= 0.0001 || 1 - glm::dot(newNormal, hitInfo.normal) * glm::dot(newNormal, hitInfo.normal) <= 0.0001) {
            realised--;
            continue;
        }
        auto newRay = Ray(refletedRay.origin, newNormal);
        finalColor += getFinalColor(scene, bvh, newRay, features, depth + 1);
    }

    if (DRAW_GLOSSY_CONE) {
        auto cornerA = glm::normalize(refletedRay.direction + (-a / 2) * U + (-a / 2) * V);
        auto cornerB = glm::normalize(refletedRay.direction + (-a / 2) * U + (-a / 2 + a) * V);
        auto cornerC = glm::normalize(refletedRay.direction + (-a / 2 + a) * U + (-a / 2) * V);
        auto cornerD = glm::normalize(refletedRay.direction + (-a / 2 + a) * U + (-a / 2 + a) * V);
        drawRay(Ray(refletedRay.origin, cornerA, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin, cornerB, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin, cornerC, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin, cornerD, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin + cornerA * 0.5f, refletedRay.origin + cornerB - refletedRay.origin - cornerA, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin + cornerA * 0.5f, refletedRay.origin + cornerC - refletedRay.origin - cornerA, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin + cornerB * 0.5f, refletedRay.origin + cornerD - refletedRay.origin - cornerB, 0.5), glm::vec3(1, 1, 0));
        drawRay(Ray(refletedRay.origin + cornerC * 0.5f, refletedRay.origin + cornerD - refletedRay.origin - cornerC, 0.5), glm::vec3(1, 1, 0));
    }
    return finalColor / glm::vec3(realised);
}