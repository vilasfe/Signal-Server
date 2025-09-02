#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "image-ppm.hh"
#include "image.hh"

ImagePPM::ImagePPM(const size_t width, const size_t height)
  : Image(width, height, IMAGE_RGB, IMAGE_PPM)
{
  set_extension(".ppm");

  const size_t buf_size = width * height * RGB_SIZE;

  /* Allocate the canvas buffer */
  resize_canvas(buf_size);
}

auto ImagePPM::add_pixel(const uint8_t r,const uint8_t g,const uint8_t b) -> int {
	set_canvas_item(get_next_pixel(), r);
	set_canvas_item(get_next_pixel()+1, g);
	set_canvas_item(get_next_pixel()+2, b);

	set_next_pixel_index(get_next_pixel() + 3);

	return 0;
}

auto ImagePPM::add_pixel(const uint8_t r, const uint8_t g, const uint8_t b, [[maybe_unused]] const uint8_t a) -> int {
	return add_pixel(r, g, b);
}
