#include "texture.h"
#include <framework/image.h>
#include <glm/common.hpp>

float clamp(float value, float min, float max) {
    return glm::min(max, glm::max(min, value));
}

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    float pixelsC = image.width * image.height - 1;
    if (features.extra.enableBilinearTextureFiltering) {
        float maxCoord = image.width * image.height;
        float u_p = fmod(texCoord.x, 1.0) * image.width - 0.5;
        float v_p = fmod(texCoord.y, 1.0) * image.height - 0.5;
        float iu0 = std::max(0.0f, floor(u_p));
        float iu1 = std::min(maxCoord, iu0 + 1);
        float iv0 = std::max(0.0f, image.width - floor(v_p));
        float iv1 = std::min(maxCoord, iv0 - 1);
        float a_u = (iu1 - u_p);
        float b_u = 1 - a_u;
        float a_v = (floor(v_p) + 1 - v_p);
        float b_v = 1 - a_v;
        return a_u * a_v * image.pixels[clamp(fmod(iu0 + image.width * iv0, maxCoord),0,pixelsC)] + a_u * b_v * image.pixels[clamp(fmod(iu0 + image.width * iv1, maxCoord),0,pixelsC)] + b_u * a_v * image.pixels[clamp(fmod(iu1 + image.width * iv0, maxCoord),0,pixelsC)] + b_u * b_v * image.pixels[clamp(fmod(iu1 + image.width * iv1, maxCoord),0,pixelsC)];
    }
    
    int w = round(texCoord.x * image.width - 0.5);
    int h = image.height - round(texCoord.y * image.height + 0.5);

    return image.pixels[clamp(h * image.width + w, 0,pixelsC)];
    
    
}