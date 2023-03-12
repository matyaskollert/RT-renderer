#pragma once
#include <framework/ray.h>
#include "../intersect.h"
#include <glm/geometric.hpp>
#include <iostream>
#include "../interpolate.h"
extern float FORCE_MIPMAP_LEVEL;
extern bool FORCE_MIPMAP;
extern int LARGEST_TEXTURE_LEVEL;

float findMipMapLevel(HitInfo hitInfo, Ray ray, const Features& features);
std::vector<std::shared_ptr<Image>> applyMipMap(std::shared_ptr<Image> img);