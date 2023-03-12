#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include "draw.h"
#include "extra/bloom_effect.h"
#include "extra/depth_of_field.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include "extra/mipmap.h"
#include <iostream>
#include <extra/glossy_reflection.h>

int MAX_RAYS_DEPTH = 5;
bool SHOW_ONE_RAY = false;
int SHOW_RAY_NUMBER = 0;
int TIME_RANGE = 1500;

void drawRaysForMotionBlur(Ray ray, const Features& features)
{

    for (int time = 0; time < 20; time++) {

        // Generate a ray with the origin at cameraPos, going through the given pixel
        // for each iteration, we trace a ray with different time stamps
        float alpha = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        glm::vec3 mov = { features.extra.motionX, features.extra.motionY, features.extra.motionZ };
        glm::vec finalMov = glm::normalize(mov);

        ray.origin = ray.origin + finalMov * (alpha - 0.5f) * features.extra.motionTime;
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
    }
}

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    glm::vec3 Lo(0.0f);

    // Checks that the show only one ray flag is set and that
    // this is the correct depth of ray to show.
    if (enableDebugDraw && SHOW_ONE_RAY) {
        HideDrawings = SHOW_RAY_NUMBER != rayDepth;
    }

    if (enableDebugDraw && features.extra.enableMotionBlur) {
        drawRaysForMotionBlur(ray, features);
    }
	
	
	
    if (bvh.intersect(ray, hitInfo, features)) {
        Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        // Checks that recursive tracing is enabled and that we haven't traced rays too deep.
        if (features.enableRecursive) {

            if (rayDepth >= MAX_RAYS_DEPTH) {
                if (SHOW_ONE_RAY)
                    HideDrawings = false;
            } else {
                // Checks that the ks is not black.
                if (glm::length(hitInfo.material.ks) > 1e-7) {
                    Ray reflection = computeReflectionRay(ray, hitInfo);
                    glm::vec3 reflectedColor;
                    if (features.extra.enableGlossyReflection) {
                        reflectedColor = getGlossy(scene, bvh, features, reflection, ray, hitInfo, rayDepth + 1);
                    } else {
                        reflectedColor = getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
                    }
                    Lo = Lo + hitInfo.material.ks * reflectedColor;

                    if (enableDebugDraw && SHOW_ONE_RAY) {
                        HideDrawings = SHOW_RAY_NUMBER != rayDepth;
                    }
                }
            }
        }

        // Transparency extra feature
        if (features.extra.enableTransparency && hitInfo.material.transparency < (1 - 0.00001f)) {
            if (rayDepth >= MAX_RAYS_DEPTH) {
                if (SHOW_ONE_RAY)
                    HideDrawings = false;
            } else {
                Ray continuation(ray.origin + ray.direction * (ray.t + 0.00001f), ray.direction);
                auto transparentColor = getFinalColor(scene, bvh, continuation, features, rayDepth + 1);
                Lo = (hitInfo.material.transparency) * Lo + (1 - hitInfo.material.transparency) * transparentColor;

                if (enableDebugDraw && SHOW_ONE_RAY) {
                    HideDrawings = SHOW_RAY_NUMBER != rayDepth;
                }
            }
        }
        // Draw a white debug ray if the ray hits.
        if (features.enableShading && DrawShadingRay) {
            drawRay(ray, Lo);
        } else {
            drawRay(ray, glm::vec3(1));
        }
    } else {
        if (features.extra.enableEnvironmentMapping) {
            glm::vec3 texture = environment_mapping(ray.direction);
            std::shared_ptr<Image> image;
            glm::vec2 texCoord = { texture.y, texture.z };
            texCoord.x = glm::clamp(texCoord.x, 0.0f, 1.0f);
            texCoord.y = glm::clamp(texCoord.y, 0.0f, 1.0f);
            if (texture.x == 0) {
                image = scene.environmentMapLeft;
            } else if (texture.x == 1) {
                image = scene.environmentMapRight;
            } else if (texture.x == 2) {
                image = scene.environmentMapTop;
            } else if (texture.x == 3) {
                image = scene.environmentMapBottom;
            } else if (texture.x == 4) {
                image = scene.environmentMapFront;
            } else if (texture.x == 5) {
                image = scene.environmentMapBack;
            }
            drawRay(ray, Lo);
            Lo = acquireTexel(*image.get(), texCoord, features);
        } else {
            // Draw a red debug ray if the ray missed.
            drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        }

    }

    if (features.extra.enableDepthOfField && rayDepth == 0) {
        addDepthOfField(scene, bvh, ray, features, Lo);
    }

    if (SHOW_ONE_RAY && rayDepth >= MAX_RAYS_DEPTH) {
        HideDrawings = false;
    }

    return Lo;
}


glm::vec3 multipleRays(const Screen& screen, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, glm::vec2 normalizedPixelPos, const Features& features)
{
    glm::vec3 color;
    glm::ivec2 windowResolution = screen.resolution();
    // number of rays casted for a single pixel
    const int numberOfRays = 4;

    // used for averaging the color
    float colorAvgx = 0.0f;
    float colorAvgy = 0.0f;
    float colorAvgz = 0.0f;
    glm::vec2 samplePoint;
    for (int i = 0; i < numberOfRays; i++) {
        for (int j = 0; j < numberOfRays; j++) {
            
            float alpha = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) ;
            float beta = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) ;
            // Generate a ray with the origin at cameraPos, going through the given pixel
            // for each iteration, we trace a ray for the same pixel, but with an offset
            // we want the number i and j to be small such that it will not exceed the size of one pixel
            samplePoint = { ((i + alpha)/ 4.0f) + normalizedPixelPos.x, ((j + beta) / 4.0f) + normalizedPixelPos.y};
            
			samplePoint.x = samplePoint.x / (windowResolution.x) * 2.0f - 1.0f;
            samplePoint.y = samplePoint.y / (windowResolution.y) * 2.0f - 1.0f;
			
            const Ray cameraRay = camera.generateRay(samplePoint);

            color = getFinalColor(scene, bvh, cameraRay, features, 0);
            colorAvgx += color.x;
            colorAvgy += color.y;
            colorAvgz += color.z;
        }
    }

    colorAvgx /= 16.0f;
    colorAvgy /= 16.0f;
    colorAvgz /= 16.0f;

    return glm::vec3 { colorAvgx, colorAvgy, colorAvgz };
}

// returns color of the pixel
glm::vec3 motionBlur(const Scene& scene, const Trackball& camera, const Features& features, const BvhInterface& bvh, glm::vec2 normalizedPixelPos)
{
    glm::vec3 color;

    // used for averaging the color
    float averageX = 0.0f;
    float averageY = 0.0f;
    float averageZ = 0.0f;
    for (int time = 0; time < TIME_RANGE; time++) {

        // Generate a ray with the origin at cameraPos, going through the given pixel
        // for each iteration, we trace a ray with different time stamps
        Ray cameraRay = camera.generateRay(normalizedPixelPos);
        float alpha = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        glm::vec3 mov = { features.extra.motionX, features.extra.motionY, features.extra.motionZ };
        glm::vec finalMov = glm::normalize(mov);

        cameraRay.origin = cameraRay.origin + finalMov * features.extra.motionTime * (alpha - 0.5f);
        color = getFinalColor(scene, bvh, cameraRay, features, 0);
        averageX += color.x;
        averageY += color.y;
        averageZ += color.z;
    }
    return { averageX / TIME_RANGE, averageY / TIME_RANGE, averageZ / TIME_RANGE };
}


glm::vec3 environment_mapping(glm::vec3 ray)
{
    int index;
    float u;
    float v;

	// use absolute values to see in which side would be intersected by the ray
    float absX = std::abs(ray.x);
    float absY = std::abs(ray.y);
    float absZ = std::abs(ray.z);

	// we use this to check on which side of the side out of the possible side pairs we are on 
    int isXPositive = ray.x > 0 ? 1 : 0;
    int isYPositive = ray.y > 0 ? 1 : 0;
    int isZPositive = ray.z > 0 ? 1 : 0;

    float maxAxis;
    float uc;
    float vc;
	
    if (absX >= absY && absX >= absZ) {
        // x is positive
        if (isXPositive) {
            // u (0 to 1) goes from +z to -z
            // v (0 to 1) goes from -y to +y
            maxAxis = absX;
            uc = -ray.z;
            vc = ray.y;
            index = 0;
            // Convert range from -1 to 1 to 0 to 1
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        } 
        // x is negative
        else {
            // u (0 to 1) goes from -z to +z
            // v (0 to 1) goes from -y to +y
            maxAxis = absX;
            uc = ray.z;
            vc = ray.y;
            index = 1;
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        }
    }
    
    if (absY >= absX && absY >= absZ) {
		
		// y is positive
        if (isYPositive) {
            // u (0 to 1) goes from -x to +x
            // v (0 to 1) goes from +z to -z
            maxAxis = absY;
            uc = ray.x;
            vc = -ray.z;
            index = 2;
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        }
		// y is negative
		else {
            // u (0 to 1) goes from -x to +x
            // v (0 to 1) goes from -z to +z
            maxAxis = absY;
            uc = ray.x;
            vc = ray.z;
            index = 3;
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        }
        
    }
    
    if (absZ >= absX && absZ >= absY) {
		// positive z
        if (isZPositive) {
            // u (0 to 1) goes from -x to +x
            // v (0 to 1) goes from -y to +y
            maxAxis = absZ;
            uc = ray.x;
            vc = ray.y;
            index = 4;
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        } 
        // negative z
        else {
            // u (0 to 1) goes from +x to -x
            // v (0 to 1) goes from -y to +y
            maxAxis = absZ;
            uc = -ray.x;
            vc = ray.y;
            index = 5;
            u = 0.5f * (uc / maxAxis + 1.0f);
            v = 0.5f * (vc / maxAxis + 1.0f);
            return { index, u, v };
        }
        
    }
}



void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            glm::vec2 pixel = normalizedPixelPos;
			 // if motion blur is activated
            if (features.extra.enableMotionBlur == true) {
                screen.setPixel(x, y, motionBlur(scene, camera, features, bvh, pixel));
            } else {
                if (features.extra.enableMultipleRaysPerPixel) {
                    glm::vec2 pix1 = { float(x), float(y) };
                    glm::vec3 colorMultipleRays = multipleRays(screen, scene, camera, bvh, pix1, features);
                    screen.setPixel(x, y, colorMultipleRays);
                } else {
                    // else compute the colors normally
                    const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                    screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features, 0));
                }
            }
        }
    }

    if (features.extra.enableBloomEffect) {
        applyBloomEffect(screen);
    }

}
