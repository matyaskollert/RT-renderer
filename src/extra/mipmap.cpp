#include <framework/ray.h>
#include "../intersect.h"
#include <glm/geometric.hpp>
#include <iostream>
#include "../interpolate.h"
#include "mipmap.h"
#include <draw.h>
bool FORCE_MIPMAP = false;
float FORCE_MIPMAP_LEVEL = 0;
int LARGEST_TEXTURE_LEVEL = 0;
glm::vec3 applyBoxFilter(std::shared_ptr<Image> img, int x, int y)
{
    glm::vec3 new_color = 
        img->pixels[y * img->width + x] + 
        img->pixels[(y+1) * img->width + x] + 
        img->pixels[y * img->width + x + 1] + 
        img->pixels[(y+1) * img->width + x + 1];
    return new_color * glm::vec3(1/4.0f);
}

std::shared_ptr<Image> applyHalving(std::shared_ptr<Image> img)
{
    std::vector<glm::vec3> pixels;
    int width = img->width / 2;
    int height = img->height / 2;

    for (int y = 0; y < img->height; y += 2) {
        for (int x = 0; x < img->width; x += 2) {
            if (x + 1 >= img->width || y + 1 >= img->height) {
                continue;
            }
            auto coord = y * img->width + x;
            glm::vec3 pixel = img->pixels[coord];
            auto averaged = applyBoxFilter(img, x, y);
            pixels.push_back(averaged);
        }
    }
    return std::make_shared<Image> (Image(width, height, pixels));
}

std::vector<std::shared_ptr<Image>> applyMipMap(std::shared_ptr<Image> img)
{
    std::vector<std::shared_ptr<Image>> results;
    int width = img->width;
    int height = img->height;
    int levels = (int)glm::log2((float)width);
    LARGEST_TEXTURE_LEVEL = glm::max(LARGEST_TEXTURE_LEVEL, levels);
    results.push_back(img);
    while (levels--) {
        results.push_back(applyHalving(results[results.size()-1]));
    }

    return results;
}

float findMipMapLevel(HitInfo hitInfo, Ray ray, const Features& features)
{
    int maxLevel = hitInfo.material.mipmaps.size() - 1;
    if (features.extra.enableMipmapTextureFiltering) {
        if (FORCE_MIPMAP)
            return glm::min(FORCE_MIPMAP_LEVEL, (float)maxLevel);

        return glm::min(glm::max(0.0f,glm::log2(ray.t)), (float)maxLevel);
    } else {
        return 0;
    }
}