#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <extra/mipmap.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 finalColor = glm::vec3(0);
    finalColor += phongSpecularOnly(hitInfo, lightPosition, (ray.origin + ray.direction * ray.t), ray.origin, lightColor);
    finalColor += diffuseOnly(hitInfo, lightPosition, (ray.origin + ray.direction * ray.t), lightColor, features, ray);
    return finalColor;
}

glm::vec3 diffuseOnly(const HitInfo hitInfo, const glm::vec3& lPos, const glm::vec3& vPos, const glm::vec3& lightColor, const Features& features, Ray ray)
{
    glm::vec3 lightPos = lPos - vPos;
    float dot = glm::dot(glm::normalize(lightPos), hitInfo.normal);
    if (dot < 0) {
        return glm::vec3(0);
    }
    if (features.enableTextureMapping && hitInfo.material.kdTexture) {
        glm::vec3 c;
        if (features.extra.enableMipmapTextureFiltering) {
            float mipmapLevel = findMipMapLevel(hitInfo, ray, features);

            // Doing trilinear interpolation
            auto c1 = acquireTexel(*hitInfo.material.mipmaps[floor(mipmapLevel)], hitInfo.texCoord, features);
            auto secondLevel = glm::min((int)hitInfo.material.mipmaps.size() - 1, (int)floor(mipmapLevel) + 1);
            auto c2 = acquireTexel(*hitInfo.material.mipmaps[secondLevel], hitInfo.texCoord, features);
            c = c1 * (floor(mipmapLevel) + 1 - mipmapLevel) + c2 * (mipmapLevel - floor(mipmapLevel));
        } else {
            c = acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
        }
        glm::vec3 mul { lightColor.x * c.x, lightColor.y * c.y, lightColor.z * c.z };
        return mul * dot;
    }
    glm::vec3 mul { lightColor.x * hitInfo.material.kd.x, lightColor.y * hitInfo.material.kd.y, lightColor.z * hitInfo.material.kd.z };
    return mul * dot;
}

glm::vec3 phongSpecularOnly(const HitInfo hitInfo, const glm::vec3& lPos, const glm::vec3& vPos, const glm::vec3& cPos, const glm::vec3& lightColor)
{
    glm::vec3 N = hitInfo.normal;
    glm::vec3 L = vPos - lPos;
    glm::vec R = L - 2 * glm::dot(L, N) * N;
    float dot = glm::dot(glm::normalize(R), glm::normalize(cPos - vPos));
    if (dot <= 0 || glm::dot(L, N) >= 0) {
        return glm::vec3(0);
    }
    glm::vec3 mul { lightColor.x * hitInfo.material.ks.x, lightColor.y * hitInfo.material.ks.y, lightColor.z * hitInfo.material.ks.z };
    return mul * glm::pow(dot, hitInfo.material.shininess);
}

const Ray computeReflectionRay(Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    auto d = glm::normalize(ray.direction);
    auto reflectionDirection = d - 2 * glm::dot(hitInfo.normal, d) * hitInfo.normal;
    Ray reflectionRay(ray.origin + ray.direction * ray.t, reflectionDirection);

    return reflectionRay;
}