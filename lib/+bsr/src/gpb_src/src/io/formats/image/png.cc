/*
 * Png file utilities.
 */
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "io/formats/image/png.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_file_open_error.hh"
#include "lang/exceptions/ex_file_read_error.hh"
#include "lang/exceptions/ex_file_write_error.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/string.hh"
#include "lang/types/type_ranges.hh"
#include "lang/types/type_sizes.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"

#include <stdio.h>

extern "C" {
#include <png.h>
} /* extern "C" */

namespace io {
namespace formats {
namespace image {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::exceptions::exception;
using lang::exceptions::ex_file_open_error;
using lang::exceptions::ex_file_read_error;
using lang::exceptions::ex_file_write_error;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;
using lang::string;
using math::matrices::matrix;

namespace {
/*
 * Open the specified file and initialize png data structures for reading.
 * Return a pointer to the file.
 */
FILE* png_initialize_read(
   const string<>& filename,
   png_structp&    png_ptr,
   png_infop&      info_ptr)
{
   /* open file for reading */
   FILE* f = fopen(filename, "r");
   if (f == NULL)
      throw ex_file_open_error(filename, "read");
   /* initialize png data structure */
   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (png_ptr == NULL) {
      fclose(f);
      throw exception("unable to create png_struct data structure");
   }
   /* initialize png info data structure */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      fclose(f);
      throw exception("unable to create png_info data structure");
   }
   /* setup input control */
   png_init_io(png_ptr, f);
   /* return file pointer */
   return f;
}

/*
 * Open the specified file and initialize png data structures for writing.
 * Return a pointer to the file.
 */
FILE* png_initialize_write(
   const string<>& filename,
   png_structp&    png_ptr,
   png_infop&      info_ptr)
{
   /* open file for writing */
   FILE* f = fopen(filename, "w");
   if (f == NULL)
      throw ex_file_open_error(filename, "write");
   /* initialize png data structure */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (png_ptr == NULL) {
      fclose(f);
      throw exception("unable to create png_struct data structure");
   }
   /* initialize png info data structure */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      png_destroy_write_struct(&png_ptr, NULL);
      fclose(f);
      throw exception("unable to create png_info data structure");
   }
   /* setup output control */
   png_init_io(png_ptr, f);
   /* return file pointer */
   return f;
}

/*
 * Allocate memory for png image data.
 */
png_bytepp create_row_pointers(png_uint_32 height, png_uint_32 row_bytes) {
   png_bytepp row_pointers = new png_bytep[height];
   for (unsigned long n = 0; n < static_cast<unsigned long>(height); n++)
      row_pointers[n] = new png_byte[row_bytes];
   return row_pointers;
}

/*
 * Deallocate memory for png image data.
 */
void destroy_row_pointers(png_bytepp row_pointers, png_uint_32 height) {
   for (unsigned long n = 0; n < static_cast<unsigned long>(height); n++)
      delete [] (row_pointers[n]);
   delete [] row_pointers;
}

/*
 * Extract real-valued image channels from a png data array.
 * Return one matrix per channel.
 */
auto_collection< matrix<>, array_list< matrix<> > > decode_image(
   png_bytepp    row_pointers,
   png_uint_32   height, 
   png_uint_32   width, 
   int           bit_depth, 
   unsigned long n_channels)
{
   /* allocate matrix for each channel */
   auto_collection< matrix<>, array_list< matrix<> > > channels(
      new array_list< matrix<> >()
   );
   for (unsigned long c = 0; c < n_channels; c++) {
      auto_ptr< matrix<> > data(new matrix<>(height, width));
      channels->add(*data);
      data.release();
   }
   /* compute step size between array entries */
   unsigned long step = static_cast<unsigned long>(bit_depth) / BITS_PER_BYTE;
   if ((step * BITS_PER_BYTE) < static_cast<unsigned long>(bit_depth))
      step++;
   unsigned long row_bytes =
      static_cast<unsigned long>(width) * step * n_channels;
   /* extract image data */
   for (unsigned long n = 0, x = 0; x < static_cast<unsigned long>(height); x++)
   {
      png_bytep row = row_pointers[x];
      for (unsigned long y = 0; y < row_bytes; n++) {
         for (unsigned long c = 0; c < n_channels; c++, y += step) {
            if (bit_depth == 8) {
               unsigned char* val = reinterpret_cast<unsigned char*>(&row[y]);
               ((*channels)[c])[n] = static_cast<double>(*val) / UCHAR_MAX;
            } else if (bit_depth == 16) {
               unsigned short* val = reinterpret_cast<unsigned short*>(&row[y]);
               ((*channels)[c])[n] = static_cast<double>(*val) / USHRT_MAX;
            } else {
               throw ex_invalid_argument("invalid bit depth specified");
            }
         }
      }
   }
   return channels;
}

/*
 * Encode real-valued image channels into a png data array.
 * Return the png image data.
 */
png_bytepp encode_image(
   const array_list< const matrix<> > channels,
   unsigned long                      bit_depth)
{
   /* get number of channels */
   unsigned long n_channels = channels.size();
   if (n_channels == 0)
      return NULL;
   /* check that all channels are two dimensional */
   for (unsigned long c = 0; c < n_channels; c++) {
      if (channels[c].dimensionality() != 2)
         throw ex_invalid_argument("image channels must be 2D");
   }
   /* get image dimensions */
   unsigned long height = channels[0].size(0);
   unsigned long width  = channels[0].size(1);
   /* check that all channels have the same dimensions */
   for (unsigned long c = 1; c < n_channels; c++) {
      if ((channels[c].size(0) != height) || (channels[c].size(1) != width))
         throw ex_invalid_argument("image channels must have equal dimensions");
   }
   /* compute step size between array entries */
   unsigned long step = bit_depth / BITS_PER_BYTE;
   if ((step * BITS_PER_BYTE) < bit_depth) { step++; }
   unsigned long row_bytes = width * step * n_channels;
   /* allocate png image */
   png_bytepp row_pointers = create_row_pointers(
      static_cast<png_uint_32>(height), static_cast<png_uint_32>(row_bytes)
   );
   /* encode image data */
   for (unsigned long n = 0, x = 0; x < height; x++) {
      png_bytep row = row_pointers[x];
      for (unsigned long y = 0; y < row_bytes; n++) {
         for (unsigned long c = 0; c < n_channels; c++, y += step) {
            if (bit_depth == 8) {
               double val = math::round((channels[c])[n] * UCHAR_MAX);
               unsigned char* im = reinterpret_cast<unsigned char*>(&row[y]);
               *im = static_cast<unsigned char>(val);
            } else if (bit_depth == 16) {
               double val = math::round((channels[c])[n] * USHRT_MAX);
               unsigned short* im = reinterpret_cast<unsigned short*>(&row[y]);
               *im = static_cast<unsigned short>(val);
            } else {
               destroy_row_pointers(
                  row_pointers, static_cast<png_uint_32>(height)
               );
               throw ex_invalid_argument("invalid bit depth specified");
            }
         }
      }
   }
   return row_pointers;
}
} /* namespace */

/*
 * Check whether the given png file contains a grayscale or color image.
 */
bool png::is_grayscale(const string<>& filename) {
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* cleanup */
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return whether image is grayscale */
   return ((color_type == PNG_COLOR_TYPE_GRAY) ||
           (color_type == PNG_COLOR_TYPE_GRAY_ALPHA));
}

/*
 * Check whether the given png file contains an alpha (transparency) channel.
 */
bool png::has_alpha_channel(const string<>& filename) {
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* cleanup */
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return whether image is grayscale */
   return ((color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ||
           (color_type == PNG_COLOR_TYPE_RGB_ALPHA));
}

/*
 * Read a grayscale image from a png file (discard any alpha channel).
 * The returned image has values in [0,1].
 * Throw an exception if the file contains a color image.
 */
void png::read(
   const string<>&       filename,
   auto_ptr< matrix<> >& image_gray)
{
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* set options for extracting grayscale image */
   if (color_type == PNG_COLOR_TYPE_GRAY) {
      /* expand grayscale image to at least 8 bits per pixel */
      if (bit_depth < 8) {
         png_set_expand_gray_1_2_4_to_8(png_ptr);
         png_read_update_info(png_ptr, info_ptr);
         bit_depth = 8;
      }
   } else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
      /* strip alpha channel */
      png_set_strip_alpha(png_ptr);
      png_read_update_info(png_ptr, info_ptr);
   } else {
      /* error - file doesn't contain grayscale image */
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(f);
      throw ex_file_read_error(filename, "grayscale image not found");
   }
   /* allocate memory for image data */
   png_bytepp row_pointers = create_row_pointers(
      height, png_get_rowbytes(png_ptr, info_ptr)
   );
   /* read image data */
   png_read_image(png_ptr, row_pointers);
   /* extract image into real-valued matrix */
   auto_collection< matrix<>, array_list< matrix<> > > im = decode_image(
      row_pointers, height, width, bit_depth, 1
   );
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return image */
   image_gray.reset(&(im->remove_head()));
}

/*
 * Read a grayscale image with alpha channel from a png file.
 * The returned grayscale and alpha channels have values in [0,1].
 * Throw an exception if the file contains a color image.
 *
 * Note that if the file does not define an alpha channel, the
 * returned alpha channel is set to one (completely opaque).
 */
void png::read(
   const string<>&       filename,
   auto_ptr< matrix<> >& image_gray,
   auto_ptr< matrix<> >& image_alpha)
{
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* set options for extracting grayscale image */
   if (color_type == PNG_COLOR_TYPE_GRAY) {
      /* expand grayscale image to at least 8 bits per pixel */
      if (bit_depth < 8) {
         png_set_expand_gray_1_2_4_to_8(png_ptr);
         png_read_update_info(png_ptr, info_ptr);
         bit_depth = 8;
      }
   } else if (color_type != PNG_COLOR_TYPE_GRAY_ALPHA) {
      /* error - file doesn't contain grayscale image */
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(f);
      throw ex_file_read_error(filename, "grayscale image not found");
   }
   /* allocate memory for image data */
   png_bytepp row_pointers = create_row_pointers(
      height, png_get_rowbytes(png_ptr, info_ptr)
   );
   /* read image data */
   png_read_image(png_ptr, row_pointers);
   /* extract image channels into real-valued matrices */
   bool has_alpha = (color_type == PNG_COLOR_TYPE_GRAY_ALPHA);
   auto_collection< matrix<>, array_list< matrix<> > > im = decode_image(
      row_pointers, height, width, bit_depth, (has_alpha ? 2 : 1)
   );
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return image */
   image_gray.reset(&(im->remove_head()));
   if (has_alpha) {
      image_alpha.reset(&(im->remove_head()));
   } else {
      image_alpha.reset(
         new matrix<>(
            static_cast<unsigned long>(height),
            static_cast<unsigned long>(width), 1.0
         )
      );
   }
}

/*
 * Read an RGB color image from a png file (discard any alpha channel).
 * The returned red, green, and blue channels have values in [0,1].
 * Throw an exception if the file contains a grayscale image.
 */
void png::read(
   const string<>&       filename,
   auto_ptr< matrix<> >& image_r,
   auto_ptr< matrix<> >& image_g,
   auto_ptr< matrix<> >& image_b)
{
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* set options for extracting rgb image */
   if (color_type == PNG_COLOR_TYPE_PALETTE) {
      /* convert paletted image to rgb */
      png_set_palette_to_rgb(png_ptr);
      png_read_update_info(png_ptr, info_ptr);
      bit_depth = 8;
   } else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) { 
      /* strip alpha channel */
      png_set_strip_alpha(png_ptr);
      png_read_update_info(png_ptr, info_ptr);
   } else if (color_type != PNG_COLOR_TYPE_RGB) {
      /* error - file doesn't contain color image */
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(f);
      throw ex_file_read_error(filename, "color image not found");
   }
   /* allocate memory for image data */
   png_bytepp row_pointers = create_row_pointers(
      height, png_get_rowbytes(png_ptr, info_ptr)
   );
   /* read image data */
   png_read_image(png_ptr, row_pointers);
   /* extract image channels into real-valued matrices */
   auto_collection< matrix<>, array_list< matrix<> > > im = decode_image(
      row_pointers, height, width, bit_depth, 3
   );
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return image */
   image_r.reset(&(im->remove_head()));
   image_g.reset(&(im->remove_head()));
   image_b.reset(&(im->remove_head()));
}

/*
 * Read an RGBA (RGB with alpha) color image from a png file.
 * The returned red, green, blue, and alpha channels have values in [0,1].
 * Throw an exception if the file contains a grayscale image.
 *
 * Note that if the file does not define an alpha channel, the
 * returned alpha channel is set to one (completely opaque).
 */
void png::read(
   const string<>&       filename,
   auto_ptr< matrix<> >& image_r,
   auto_ptr< matrix<> >& image_g,
   auto_ptr< matrix<> >& image_b,
   auto_ptr< matrix<> >& image_a)
{
   /* open file and initialize png data for read */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f = png_initialize_read(filename, png_ptr, info_ptr);
   /* read file info */
   png_read_info(png_ptr, info_ptr);
   /* extract fields from file info */
   png_uint_32 width, height;
   int bit_depth, color_type, interlace_method, compress_method, filter_method;
   png_get_IHDR(
      png_ptr, info_ptr,
      &width, &height,
      &bit_depth, &color_type,
      &interlace_method, &compress_method, &filter_method
   );
   /* set options for extracting rgb image */
   if (color_type == PNG_COLOR_TYPE_PALETTE) {
      /* convert paletted image to rgb */
      png_set_palette_to_rgb(png_ptr);
      png_read_update_info(png_ptr, info_ptr);
      bit_depth = 8;
   } else if ((color_type != PNG_COLOR_TYPE_RGB) &&
              (color_type != PNG_COLOR_TYPE_RGB_ALPHA)) 
   {
      /* error - file doesn't contain color image */
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(f);
      throw ex_file_read_error(filename, "color image not found");
   }
   /* allocate memory for image data */
   png_bytepp row_pointers = create_row_pointers(
      height, png_get_rowbytes(png_ptr, info_ptr)
   );
   /* read image data */
   png_read_image(png_ptr, row_pointers);
   /* extract image channels into real-valued matrices */
   bool has_alpha = (color_type == PNG_COLOR_TYPE_RGB_ALPHA);
   auto_collection< matrix<>, array_list< matrix<> > > im = decode_image(
      row_pointers, height, width, bit_depth, (has_alpha ? 4 : 3)
   );
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
   fclose(f);
   /* return image */
   image_r.reset(&(im->remove_head()));
   image_g.reset(&(im->remove_head()));
   image_b.reset(&(im->remove_head()));
   if (has_alpha) {
      image_a.reset(&(im->remove_head()));
   } else {
      image_a.reset(
         new matrix<>(
            static_cast<unsigned long>(height),
            static_cast<unsigned long>(width), 1.0
         )
      );
   }
}

/*
 * Write a grayscale image to a png file.
 * The image should have values in [0,1].
 */
void png::write(
   const string<>& filename,
   const matrix<>& image_gray,
   unsigned long   bit_depth)
{
   /* encode image */
   array_list< const matrix<> > channels;
   channels.add(image_gray);
   png_bytepp row_pointers = encode_image(channels, bit_depth);
   /* get image size */
   png_uint_32 height = static_cast<png_uint_32>(image_gray.size(0));
   png_uint_32 width  = static_cast<png_uint_32>(image_gray.size(1));
   /* open file and initialize png data for write */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f;
   try {
      f = png_initialize_write(filename, png_ptr, info_ptr);
   } catch (...) {
      destroy_row_pointers(row_pointers, height);
      throw;
   }  
   /* set image information */
   png_set_IHDR(
      png_ptr, info_ptr, 
      width, height,
      static_cast<int>(bit_depth),
      PNG_COLOR_TYPE_GRAY,
      PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE,
      PNG_FILTER_TYPE_BASE
   );
   /* write image information */
   png_write_info(png_ptr, info_ptr);
   /* write image data */
   png_write_image(png_ptr, row_pointers);
   png_write_end(png_ptr, info_ptr);
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(f);
}

/*
 * Write a grayscale image with alpha channel to a png file.
 * The grayscale and alpha channels should have values in [0,1].
 */
void png::write(
   const string<>& filename,
   const matrix<>& image_gray,
   const matrix<>& image_alpha,
   unsigned long   bit_depth)
{
   /* encode image */
   array_list< const matrix<> > channels;
   channels.add(image_gray);
   channels.add(image_alpha);
   png_bytepp row_pointers = encode_image(channels, bit_depth);
   /* get image size */
   png_uint_32 height = static_cast<png_uint_32>(image_gray.size(0));
   png_uint_32 width  = static_cast<png_uint_32>(image_gray.size(1));
   /* open file and initialize png data for write */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f;
   try {
      f = png_initialize_write(filename, png_ptr, info_ptr);
   } catch (...) {
      destroy_row_pointers(row_pointers, height);
      throw;
   }  
   /* set image information */
   png_set_IHDR(
      png_ptr, info_ptr, 
      width, height,
      static_cast<int>(bit_depth),
      PNG_COLOR_TYPE_GRAY_ALPHA,
      PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE,
      PNG_FILTER_TYPE_BASE
   );
   /* write image information */
   png_write_info(png_ptr, info_ptr);
   /* write image data */
   png_write_image(png_ptr, row_pointers);
   png_write_end(png_ptr, info_ptr);
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(f);
}

/* 
 * Write an RGB color image to a png file.
 * The red, green, and blue channels should have values in [0,1].
 */
void png::write(
   const string<>& filename,
   const matrix<>& image_r,
   const matrix<>& image_g,
   const matrix<>& image_b,
   unsigned long   bit_depth)
{
   /* encode image */
   array_list< const matrix<> > channels;
   channels.add(image_r);
   channels.add(image_g);
   channels.add(image_b);
   png_bytepp row_pointers = encode_image(channels, bit_depth);
   /* get image size */
   png_uint_32 height = static_cast<png_uint_32>(image_r.size(0));
   png_uint_32 width  = static_cast<png_uint_32>(image_r.size(1));
   /* open file and initialize png data for write */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f;
   try {
      f = png_initialize_write(filename, png_ptr, info_ptr);
   } catch (...) {
      destroy_row_pointers(row_pointers, height);
      throw;
   }  
   /* set image information */
   png_set_IHDR(
      png_ptr, info_ptr, 
      width, height,
      static_cast<int>(bit_depth),
      PNG_COLOR_TYPE_RGB,
      PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE,
      PNG_FILTER_TYPE_BASE
   );
   /* write image information */
   png_write_info(png_ptr, info_ptr);
   /* write image data */
   png_write_image(png_ptr, row_pointers);
   png_write_end(png_ptr, info_ptr);
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(f);
}

/* 
 * Write an RGBA (RGB with alpha) color image to a png file.
 * The red, green, blue, and alpha channels should have values in [0,1].
 */
void png::write(
   const string<>& filename,
   const matrix<>& image_r,
   const matrix<>& image_g,
   const matrix<>& image_b,
   const matrix<>& image_a,
   unsigned long   bit_depth)
{
   /* encode image */
   array_list< const matrix<> > channels;
   channels.add(image_r);
   channels.add(image_g);
   channels.add(image_b);
   channels.add(image_a);
   png_bytepp row_pointers = encode_image(channels, bit_depth);
   /* get image size */
   png_uint_32 height = static_cast<png_uint_32>(image_r.size(0));
   png_uint_32 width  = static_cast<png_uint_32>(image_r.size(1));
   /* open file and initialize png data for write */
   png_structp png_ptr;
   png_infop info_ptr;
   FILE* f;
   try {
      f = png_initialize_write(filename, png_ptr, info_ptr);
   } catch (...) {
      destroy_row_pointers(row_pointers, height);
      throw;
   }  
   /* set image information */
   png_set_IHDR(
      png_ptr, info_ptr, 
      width, height,
      static_cast<int>(bit_depth),
      PNG_COLOR_TYPE_RGB_ALPHA,
      PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE,
      PNG_FILTER_TYPE_BASE
   );
   /* write image information */
   png_write_info(png_ptr, info_ptr);
   /* write image data */
   png_write_image(png_ptr, row_pointers);
   png_write_end(png_ptr, info_ptr);
   /* cleanup */
   destroy_row_pointers(row_pointers, height);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(f);
}

} /* namespace image */
} /* namespace formats */
} /* namespace io */
