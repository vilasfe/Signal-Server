/*
 * Generic image output handling. This feature allows
 * for extensible image output formats and permits for
 * cleaner image rendering code in the model generation
 * routines by moving all of the file-format specific
 * logic to be handled here.
 */
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring> // still using memset for some reason
#include <format>
#include <memory>
#include <print>
#include <stdexcept>
#include <string>

#include "image-ppm.hh"
#include "image.hh"

image_format Image::default_format = IMAGE_PPM;
std::string Image::dynamic_backend;

/*
 * image_set_format
 * Changes the default format for the next
 * uninitialized image canvas
 */
auto Image::set_format(image_format format) -> int {
	if(format <= IMAGE_DEFAULT || format >= IMAGE_FORMAT_MAX) {
		return EINVAL;
	}
	default_format = format;
	return 0;
}

/*
 * image_init
 * Initialize an image context. Must be called
 * before attempting to write any image data
 */
auto Image::create(const size_t width, const size_t height, [[maybe_unused]] const image_model model, const image_format format) -> std::shared_ptr<Image> {
	if(format >= IMAGE_FORMAT_MAX || (format == IMAGE_LIBRARY && dynamic_backend.empty())) {
		throw std::runtime_error("Invalid image format");
	}

	if (format == IMAGE_PPM || (format == IMAGE_DEFAULT && default_format == IMAGE_PPM)) {
		return std::make_shared<ImagePPM>(width, height);
	}
	//return std::make_shared<Image>(width, height, model, format == IMAGE_DEFAULT ? default_format : format);
	return {};
}

Image::Image(const size_t width, const size_t height, const image_model model, const image_format format)
  : width(width), height(height), model(model), format(format)

{
	/* Perform some sanity checking on provided arguments */
	if(width == 0 || height == 0) {
		throw std::runtime_error(std::format("Invalid width {} or height {} for image", width, height));
	}
	if(model < 0 || model > IMAGE_MODEL_MAX) {
		throw std::runtime_error("Invalid model parameter for image");
	}
}

void Image::set_extension(const std::string& ext) {
	extension = ext;
}

/*
 * image_add_pixel, image_set_pixel, image_get_pixel, image_write, image_free
 * Various setters ad getters for assigning pixel data. image_write
 * takes an open file handle and writes image data to it.
 * These functions simply wrap the underlying format-specific functions
 */
auto Image::write(FILE *fd) -> int {
	const size_t count = get_width() * get_height() * RGB_SIZE;

	std::println(fd, "P6\n{} {}\n255", get_width(), get_height());
	const size_t written = fwrite(canvas.data(),sizeof(uint8_t),count,fd);
	if(written < count) {
		return EPIPE;
	}

	return 0;
}

/*
 * image_get_filename
 * Creates an appropriate file name using data supplied
 * by the user. If the extension is already correct, return
 * that; if not append if there is space
 */
auto Image::get_filename(const std::string& in) const -> std::string {

	if(in.ends_with(extension)) {
		/* Already has correct extension */
		return in;
	}
	/* Doesn't have correct extension */
	return in + extension;
}

auto Image::get_width() const -> size_t {
	return width;
}

auto Image::get_height() const -> size_t {
	return height;
}

void Image::resize_canvas(size_t size) {
	canvas.clear();
	canvas.resize(size);
}

void Image::set_canvas_item(size_t index, uint8_t value) {
	canvas.at(index) = value;
}

auto Image::get_next_pixel() const -> size_t {
	return next_pixel;
}

void Image::set_next_pixel_index(size_t index) {
	next_pixel = index;
}

