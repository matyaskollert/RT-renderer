#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    glm::vec3 a = v1 - v0, b = v2 - v0, c = p - v0;
    float dotAA = glm::dot(a, a);
    float dotAB = glm::dot(a, b);
    float dotBB = glm::dot(b, b);
    float dotCA = glm::dot(c, a);
    float dotCB = glm::dot(c, b);
    float determinant = dotAA * dotBB - dotAB * dotAB;
    if (determinant == 0) {
        return false;
    }
    float deter1 = dotBB * dotCA - dotAB * dotCB;
    float alpha = deter1 / determinant;
    float deter2 = dotAA * dotCB - dotAB * dotCA;
    float beta = deter2 / determinant;
    if (alpha < 0 || beta < 0 || (alpha + beta) > 1) {
        return false;
    }
    return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float x = glm::dot(ray.direction, plane.normal);
    if (x == 0) {
        return false;
    }
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / x;
    if (t > 0) {
        if (ray.t > t) {
            ray.t = t;
        }
        return true;
    } else {
        return false;
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    plane.normal = glm::normalize(glm::cross(v0 - v2, v1 - v2));
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    float oldT = ray.t;
    if (intersectRayWithPlane(trianglePlane(v0, v1, v2), ray)) {
        if (pointInTriangle(v0, v1, v2, glm::vec3(0), (ray.origin + ray.direction * ray.t))) {

            // We need this check here. The prebuilt library does this internally.
            if (ray.t > oldT) {
                ray.t = oldT;
                return false;
            }

            return true;
        } else {
            ray.t = oldT;
        }
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    glm::vec3 o = ray.origin - sphere.center;
    float A = glm::pow(ray.direction.x, 2) + glm::pow(ray.direction.y, 2) + glm::pow(ray.direction.z, 2);
    float B = 2 * ((ray.direction.x * o.x) + (ray.direction.y * o.y) + (ray.direction.z * o.z));
    float C = glm::pow(o.x, 2) + glm::pow(o.y, 2) + glm::pow(o.z, 2) - glm::pow(sphere.radius, 2);
    float discriminant = glm::pow(B, 2) - (4 * A * C);
    float oldT = ray.t;

    if (discriminant < 0 || A == 0) {
        return false;
    } else if (discriminant == 0) {
        float t = (-B) / (2 * A);
        if (t < 0) {
            return false;
        }
        if (ray.t > t) {
            ray.t = t;
        }
    } else {
        float t1 = (-B + glm::sqrt(discriminant)) / (2 * A);
        float t2 = (-B - glm::sqrt(discriminant)) / (2 * A);
        float t;
        if (t1 < 0 && t2 < 0) {
            return false;
        } else if (t1 < 0) {
            t = t2;
        } else if (t2 < 0) {
            t = t1;
        } else {
            t = t1 < t2 ? t1 : t2;
        }
        if (ray.t > t) {
            ray.t = t;
        }
    }

    return true;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    if (ray.direction.x == 0 || ray.direction.y == 0 || ray.direction.z == 0) {
        return false;
    }
    float txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;
    float tinx = glm::min(txmin, txmax);
    float toutx = glm::max(txmin, txmax);
    float tiny = glm::min(tymin, tymax);
    float touty = glm::max(tymin, tymax);
    float tinz = glm::min(tzmin, tzmax);
    float toutz = glm::max(tzmin, tzmax);
    float tin = glm::max(tinx, glm::max(tiny, tinz));
    float tout = glm::min(toutx, glm::min(touty, toutz));
    if (tin > tout || tout < 0) {
        return false;
    }
    float t = tin < 0 ? tout : tin;
    if (ray.t > t) {
        ray.t = t;
    }
    return true;
}
