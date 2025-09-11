#ifndef IMAGE_HH_
#define IMAGE_HH_

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
  static auto create(size_t width, size_t height, _image_model model, _image_format format) -> std::shared_ptr<Image>;
  virtual ~Image() = default;
  virtual auto add_pixel(uint8_t r, uint8_t g, uint8_t b) -> int = 0;
  virtual auto add_pixel(uint8_t r, uint8_t g, uint8_t b, uint8_t a) -> int = 0;
  //virtual int get_pixel(const size_t x,const size_t y, uint8_t const* r, uint8_t const* g, uint8_t const* b, uint8_t const* a) = 0;
  [[nodiscard]] auto get_filename(const std::string& in) const -> std::string;
  virtual auto write(FILE* fd) -> int;

  Image(size_t width, size_t height, _image_model model, _image_format format);

protected:
	void set_extension(const std::string& ext);
	[[nodiscard]] auto get_width() const -> size_t;
	[[nodiscard]] auto get_height() const -> size_t;
	void resize_canvas(size_t size);
	[[nodiscard]] auto get_next_pixel() const -> size_t;
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
