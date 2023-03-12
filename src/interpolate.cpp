#include "interpolate.h"
#include <draw.h>
#include <glm/ext/matrix_float3x3.hpp>
#include <glm/geometric.hpp>
#include <iostream>

double areaOfTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    glm::vec3 a = v1 - v0;
    glm::vec3 b = v2 - v0;

    return glm::length(glm::cross(a, b));
}

glm::vec3 solve(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 P)
{
    double a1 = areaOfTriangle(P, b, c);
    double a2 = areaOfTriangle(a, P, c);
    double A = areaOfTriangle(a, b, c);

    float alpha = a1 / A;
    float beta = a2 / A;
    double gamma = 1 - alpha - beta;
    return glm::vec3 {
        alpha,
        beta,
        gamma
    };
}

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 coords = solve(v0, v1, v2, p);
    return coords;
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    auto normal = n0 * barycentricCoord.x + n1 * barycentricCoord.y + n2 * barycentricCoord.z;
    return normal;
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return t0 * barycentricCoord.x + t1 * barycentricCoord.y + t2 * barycentricCoord.z;
}