#include "screen.h"
#include <glm/gtc/random.hpp>
#include <bvh_interface.h>
#include "../draw.h"
#include "../render.h"
#include "depth_of_field.h"


float apeture = 0.1f;
float focal_length = 2.0f;
int depth_of_field_rays = 5;

// Function for adding depth of field to an image. Modifies the color in place.
void addDepthOfField(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, glm::vec3& color)
{
    HideDrawings = !DrawDepthOfField;

    // Calculates the focal point by multiplying by the focal length.
    glm::vec3 focal_point = ray.origin + ray.direction * focal_length;

    // Debug drawing of the focal point.
    drawDebugSphere(focal_point, 0.01f, glm::vec3(0, 1, 0));

    // Does x rays
    for (int i = 0; i < depth_of_field_rays; i++) {

        // Offset the origin by random number between -0.5,0.5 multiplied by apeture.
        glm::vec3 offsets = glm::vec3 {
            glm::linearRand(-0.5f * apeture, 0.5f * apeture),
            glm::linearRand(-0.5f * apeture, 0.5f * apeture),
            glm::linearRand(-0.5f * apeture, 0.5f * apeture)
        };

        // Sets new origin
        auto newOrigin = ray.origin + offsets;

        // New direction(towards the focal point from the new origin).
        auto newDirection = glm::normalize(focal_point - newOrigin);

        // Creates the new ray.
        Ray newRay = {
            newOrigin,
            newDirection,
        };

        // Adds the result of this color onto the running sum (1 as ray depth since it sufficies).
        color += getFinalColor(scene, bvh, newRay, features, 1);
    }

    // Averages the results
    color *= glm::vec3(1.0f / (depth_of_field_rays + 1.0f));
    HideDrawings = false;
}