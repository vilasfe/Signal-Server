#ifndef IMAGE_PPM_HH
#define IMAGE_PPM_HH

#include <cstdint>

#include "image.hh"

class ImagePPM : public Image {
public:
	ImagePPM(size_t width, size_t height);
	auto add_pixel(uint8_t r, uint8_t g, uint8_t b) -> int override;
	auto add_pixel(uint8_t r, uint8_t g, uint8_t b, uint8_t a) -> int override;
};

#endif
