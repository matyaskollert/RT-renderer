#include "image.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
DISABLE_WARNINGS_POP()
#include <cassert>
#include <exception>
#include <iostream>
#include <string>

Image::Image(const std::filesystem::path& filePath)
{
	if (!std::filesystem::exists(filePath)) {
		std::cerr << "Texture file " << filePath << " does not exists!" << std::endl;
		throw std::exception();
	}

	const auto filePathStr = filePath.string(); // Create l-value so c_str() is safe.
	[[maybe_unused]] int numChannelsInSourceImage;
	stbi_uc* stbPixels = stbi_load(filePathStr.c_str(), &width, &height, &numChannelsInSourceImage, STBI_rgb);

	if (!stbPixels) {
		std::cerr << "Failed to read texture " << filePath << " using stb_image.h" << std::endl;
		throw std::exception();
	}

	constexpr size_t numChannels = 3; // STBI_rgb == 3 channels
	for (size_t i = 0; i < width * height * numChannels; i += numChannels) {
            pixels.emplace_back(stbPixels[i + 0] / 255.0f, stbPixels[i + 1] / 255.0f, stbPixels[i + 2] / 255.0f);
	}

	stbi_image_free(stbPixels);
}

// Mip Map modification
Image::Image(int width, int height, std::vector<glm::vec3> pixels)
{
    this->width = width;
    this->height = height;
    this->pixels = pixels;
}