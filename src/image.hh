#ifndef _IMAGE_HH_
#define _IMAGE_HH_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#define RGB_SIZE  3
#define RGBA_SIZE 4

enum _image_format{	IMAGE_DEFAULT = 0, \
					IMAGE_PPM, \
					IMAGE_LIBRARY, \
					IMAGE_FORMAT_MAX \
				};

enum _image_model{	IMAGE_RGB, \
					IMAGE_RGBA, \
					IMAGE_MODEL_MAX
				};

class Image {
public:
  static auto set_format(_image_format format) -> int;
  static auto create(const size_t width, const size_t height, const _image_model model, const _image_format format) -> std::shared_ptr<Image>;
  virtual ~Image() = default;
  virtual int add_pixel(const uint8_t r, const uint8_t g, const uint8_t b) = 0;
  virtual int add_pixel(const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a) = 0;
  //virtual int get_pixel(const size_t x,const size_t y, uint8_t const* r, uint8_t const* g, uint8_t const* b, uint8_t const* a) = 0;
  auto get_filename(const std::string& in) const -> std::string;
  virtual int write(FILE*);

#define PIXEL_OFFSET(x,y,width,pixel_size) (((x) * (pixel_size)) + ((width) * (pixel_size) * (y)))

  Image(const size_t width, const size_t height, const _image_model model, const _image_format format);

protected:
	void set_extension(const std::string& ext);
	auto get_width() const -> size_t;
	auto get_height() const -> size_t;
	void resize_canvas(size_t size);
	auto get_next_pixel() const -> size_t;
	void set_next_pixel_index(size_t index);
	void set_canvas_item(size_t index, uint8_t value);

private:
	size_t width = 0;
	size_t height = 0;
	_image_model model = IMAGE_RGB;
	_image_format format = IMAGE_DEFAULT;
	std::vector<uint8_t> canvas;
	size_t next_pixel = 0;
	std::string extension;

	static _image_format default_format;
	static std::string dynamic_backend;
};

#endif
