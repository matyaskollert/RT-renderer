#include <screen.h>
#include "bloom_effect.h"

int FILTER_RADIUS = 1;
float BLOOM_THRESHOLD = 0.7f;
float BLOOM_SCALE = 0.2f;
bool SHOW_ONLY_BLOOM_THRESHOLD = false;
bool SHOW_ONLY_BLOOM = false;

glm::vec3 applyBoxFilter(Screen& screen, std::vector<glm::vec3>& pixels, int x, int y, int radius, float scale)
{
    glm::vec3 new_color = glm::vec3(0.0f);
    for (int y_n = y - radius; y_n <= y + radius; y_n++) {
        for (int x_n = x - radius; x_n <= x + radius; x_n++) {
            auto coord = screen.indexAt(x_n, y_n);
            if (coord < 0 || coord >= pixels.size())
                continue;

            glm::vec3& pixel = pixels[coord];
            auto percentage_of_white = (pixel.x + pixel.y + pixel.z) / 3.0;
            if (percentage_of_white >= BLOOM_THRESHOLD) {
                new_color += scale * pixel;
            }
        }
    }

    return new_color / glm::vec3((1 + 2 * radius) * (1 + 2 * radius));
}

void applyBloomEffect(Screen& screen) {
    auto pixels = screen.pixels();
    auto windowResolution = screen.resolution();


#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            auto coord = screen.indexAt(x, y);
            glm::vec3& pixel = pixels[coord];

            if (SHOW_ONLY_BLOOM_THRESHOLD) {
                auto percentage_of_white = (pixel.x + pixel.y + pixel.z) / 3.0;
                if (percentage_of_white >= BLOOM_THRESHOLD) {
                    screen.setPixel(x, y, glm::vec3(1.0));
                } else {
                    screen.setPixel(x, y, glm::vec3(0.0));
                }
            } else if (SHOW_ONLY_BLOOM) {
                auto effect = applyBoxFilter(screen, pixels, x, y, FILTER_RADIUS, BLOOM_SCALE);
                screen.setPixel(x, y, effect);
            }
            else
            {
                auto effect = applyBoxFilter(screen, pixels, x, y, FILTER_RADIUS, BLOOM_SCALE);
                screen.setPixel(x, y, pixel + effect);
            }
        }
    }
}
