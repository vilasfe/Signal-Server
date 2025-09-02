#ifndef _IMAGE_PPM_HH
#define _IMAGE_PPM_HH

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "image.hh"

class ImagePPM : public Image {
public:
	ImagePPM(const size_t width, const size_t height);
	auto add_pixel(const uint8_t r, const uint8_t g, const uint8_t b) -> int override;
	auto add_pixel(const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a) -> int override;
};

#endif
