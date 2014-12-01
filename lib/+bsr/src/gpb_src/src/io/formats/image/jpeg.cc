/*
 * Jpeg file utilities.
 */
#include "io/formats/image/jpeg.hh"
#include "lang/exceptions/ex_file_open_error.hh"
#include "lang/exceptions/ex_file_read_error.hh"
#include "lang/exceptions/ex_file_write_error.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/string.hh"
#include "lang/types/type_ranges.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"

#include <stdio.h>

extern "C" {
#include <jpeglib.h>
#include <jerror.h>
} /* extern "C" */

extern "C" {
/*
 * C wrapper for jpeg_create_compress (avoids C++ old-style cast warnings).
 */
void jpeg_init_compress(j_compress_ptr cinfo) {
   jpeg_create_compress(cinfo);
}

/*
 * C wrapper for jpeg_create_decompress (avoids C++ old-style cast warnings).
 */
void jpeg_init_decompress(j_decompress_ptr cinfo) {
   jpeg_create_decompress(cinfo);
}
} /* extern "C" */

namespace io {
namespace formats {
namespace image {
/*
 * Imports.
 */
using lang::exceptions::ex_file_open_error;
using lang::exceptions::ex_file_read_error;
using lang::exceptions::ex_file_write_error;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;
using lang::string;
using math::matrices::matrix;

/*
 * Check whether the given jpeg file contains a grayscale or color image.
 */
bool jpeg::is_grayscale(const string<>& filename) {
   /* open file for reading */
   FILE* f = fopen(filename, "r");
   if (f == NULL)
      throw ex_file_open_error(filename, "read");
   /* initialize jpeg decompression */
   struct jpeg_decompress_struct cinfo;
   struct jpeg_error_mgr jerr;
   cinfo.err = jpeg_std_error(&jerr);
   jpeg_init_decompress(&cinfo);
   jpeg_stdio_src(&cinfo, f);
   /* read jpeg file header */
   bool header_ok = (jpeg_read_header(&cinfo, true) == JPEG_HEADER_OK);
   bool is_gray = (cinfo.out_color_space == JCS_GRAYSCALE);
   jpeg_destroy_decompress(&cinfo);
   /* close file */
   fclose(f);
   /* check if header is ok */
   if (!header_ok)
      throw ex_file_read_error(filename, "could not read jpeg header");
   /* return grayscale flag */
   return is_gray;
}

/*
 * Read a grayscale image from a jpeg file.
 * The returned image has values in [0,1].
 * Throw an exception if the file contains a color image.
 */
void jpeg::read(const string<>& filename, auto_ptr< matrix<> >& image) {
   /* open file for reading */
   FILE* f = fopen(filename, "r");
   if (f == NULL)
      throw ex_file_open_error(filename, "read");
   /* initialize jpeg decompression */
   struct jpeg_decompress_struct cinfo;
   struct jpeg_error_mgr jerr;
   cinfo.err = jpeg_std_error(&jerr);
   jpeg_init_decompress(&cinfo);
   jpeg_stdio_src(&cinfo, f);
   /* check for jpeg file header */
   if (jpeg_read_header(&cinfo, true) != JPEG_HEADER_OK) {
      jpeg_destroy_decompress(&cinfo);
      fclose(f);
      throw ex_file_read_error(filename, "could not read jpeg header");
   }
   /* check for grayscale jpeg image */
   if (cinfo.out_color_space != JCS_GRAYSCALE) {
      jpeg_destroy_decompress(&cinfo);
      fclose(f);
      throw ex_file_read_error(filename, "grayscale image not found");
   }
   /* compute image dimensions */
   jpeg_calc_output_dimensions(&cinfo);
   unsigned long height = static_cast<unsigned long>(cinfo.output_height);
   unsigned long width  = static_cast<unsigned long>(cinfo.output_width);
   /* initialize image matrix and line buffer */
   auto_ptr< matrix<> > im(new matrix<>(height, width));
   unsigned char* buffer = new unsigned char[width];
   /* read image data */
   jpeg_start_decompress(&cinfo);
   for (unsigned long n = 0, x = 0; x < height; x++) {
      if (jpeg_read_scanlines(&cinfo, &buffer, 1) != 1) {
         /* error reading image data */
         jpeg_finish_decompress(&cinfo);
         jpeg_destroy_decompress(&cinfo);
         delete [] buffer;
         fclose(f);
         throw ex_file_read_error(filename, "error reading jpeg image data");
      } else {
         /* store scanline */
         for (unsigned long y = 0; y < width; y++)
            (*im)[n++] = static_cast<double>(buffer[y]) / UCHAR_MAX;
      }
   }
   /* cleanup */
   jpeg_finish_decompress(&cinfo);
   jpeg_destroy_decompress(&cinfo);
   delete [] buffer;
   fclose(f);
   /* assign output image */
   image = im;
}

/*
 * Read an RGB color image from a jpeg file.
 * The returned red, green, and blue channels have values in [0,1].
 * Throw an exception if the file contains a grayscale image.
 */
void jpeg::read(
   const string<>&       filename,
   auto_ptr< matrix<> >& image_r,
   auto_ptr< matrix<> >& image_g,
   auto_ptr< matrix<> >& image_b)
{
   /* open file for reading */
   FILE* f = fopen(filename, "r");
   if (f == NULL)
      throw ex_file_open_error(filename, "read");
   /* initialize jpeg decompression */
   struct jpeg_decompress_struct cinfo;
   struct jpeg_error_mgr jerr;
   cinfo.err = jpeg_std_error(&jerr);
   jpeg_init_decompress(&cinfo);
   jpeg_stdio_src(&cinfo, f);
   /* check for jpeg file header */
   if (jpeg_read_header(&cinfo, true) != JPEG_HEADER_OK) {
      jpeg_destroy_decompress(&cinfo);
      fclose(f);
      throw ex_file_read_error(filename, "could not read jpeg header");
   }
   /* check for color jpeg image */
   if ((cinfo.out_color_space == JCS_UNKNOWN) ||
       (cinfo.out_color_space == JCS_GRAYSCALE))
   {
      jpeg_destroy_decompress(&cinfo);
      fclose(f);
      throw ex_file_read_error(filename, "color image not found");
   }
   /* set output colorspace to RGB */
   cinfo.out_color_space = JCS_RGB;
   /* compute image dimensions */
   jpeg_calc_output_dimensions(&cinfo);
   unsigned long height = static_cast<unsigned long>(cinfo.output_height);
   unsigned long width  = static_cast<unsigned long>(cinfo.output_width);
   /* initialize image channels and line buffer */
   auto_ptr< matrix<> > im_r(new matrix<>(height, width));
   auto_ptr< matrix<> > im_g(new matrix<>(height, width));
   auto_ptr< matrix<> > im_b(new matrix<>(height, width));
   unsigned long buffer_size = 3*width; /* 3 components (RGB) per pixel */
   unsigned char* buffer = new unsigned char[buffer_size];
   /* read image data */
   jpeg_start_decompress(&cinfo);
   for (unsigned long n = 0, x = 0; x < height; x++) {
      if (jpeg_read_scanlines(&cinfo, &buffer, 1) != 1) {
         /* error reading image data */
         jpeg_finish_decompress(&cinfo);
         jpeg_destroy_decompress(&cinfo);
         delete [] buffer;
         fclose(f);
         throw ex_file_read_error(filename, "error reading jpeg image data");
      } else {
         /* store scanline */
         for (unsigned long y = 0; y < buffer_size; n++) {
            (*im_r)[n] = static_cast<double>(buffer[y++]) / UCHAR_MAX;
            (*im_g)[n] = static_cast<double>(buffer[y++]) / UCHAR_MAX;
            (*im_b)[n] = static_cast<double>(buffer[y++]) / UCHAR_MAX;
         }
      }
   }
   /* cleanup */
   jpeg_finish_decompress(&cinfo);
   jpeg_destroy_decompress(&cinfo);
   delete [] buffer;
   fclose(f);
   /* assign output color channels */
   image_r = im_r;
   image_g = im_g;
   image_b = im_b;
}

/*
 * Write a grayscale image to a jpeg file.
 * The image should have values in [0,1].
 * Optionally specify a compression quality in [0,1].
 */
void jpeg::write(
   const string<>& filename,
   const matrix<>& image,
   double          quality)
{
   /* check arguments - image dimensionality */
   if (image.dimensionality() != 2)
      throw ex_invalid_argument("image matrix must each be 2D");
   /* check arguments - quality */
   if ((quality < 0) || (quality > 1))
      throw ex_invalid_argument("image quality must be in [0,1]");
   /* open file for writing */
   FILE* f = fopen(filename, "w");
   if (f == NULL)
      throw ex_file_open_error(filename, "write");
   /* initialize jpeg compression */
   struct jpeg_compress_struct cinfo;
   struct jpeg_error_mgr jerr;
   cinfo.err = jpeg_std_error(&jerr);
   jpeg_init_compress(&cinfo);
   jpeg_stdio_dest(&cinfo, f);
   /* set image parameters */
   unsigned long height = image.size(0);
   unsigned long width  = image.size(1);
   cinfo.image_height = height;
   cinfo.image_width  = width;
   cinfo.input_components = 1;
   cinfo.in_color_space = JCS_GRAYSCALE;
   jpeg_set_defaults(&cinfo);
   jpeg_set_quality(&cinfo, static_cast<int>(quality*100), false);
   /* initialize line buffer */
   unsigned char* buffer = new unsigned char[width];
   /* write image data */
   jpeg_start_compress(&cinfo, true);
   for (unsigned long n = 0, x = 0; x < height; x++) {
      /* assemble scanline */
      for (unsigned long y = 0; y < width; y++) {
         double val = math::round(image[n++] * UCHAR_MAX);
         buffer[y] = static_cast<unsigned char>(val);
      }
      /* write scanline */
      if (jpeg_write_scanlines(&cinfo, &buffer, 1) != 1) {
         /* error writing image data */
         jpeg_finish_compress(&cinfo);
         jpeg_destroy_compress(&cinfo);
         delete [] buffer;
         fclose(f);
         throw ex_file_write_error(filename, "error writing jpeg image data");
      }
   }
   /* cleanup */
   jpeg_finish_compress(&cinfo);
   jpeg_destroy_compress(&cinfo);
   delete [] buffer;
   fclose(f);
}

/* 
 * Write an RGB color image to a jpeg file.
 * The red, green, and blue channels should have values in [0,1].
 * Optionally specify a compression quality in [0,1].
 */
void jpeg::write(
   const string<>& filename,
   const matrix<>& image_r,
   const matrix<>& image_g,
   const matrix<>& image_b,
   double          quality)
{
   /* check arguments - image dimensionality */
   if ((image_r.dimensionality() != 2) || 
       (image_g.dimensionality() != 2) || 
       (image_b.dimensionality() != 2))
      throw ex_invalid_argument("RGB image matrices must all be 2D");
   /* check arguments - image dimensions */
   unsigned long height = image_r.size(0);
   unsigned long width  = image_r.size(1);
   if ((image_g.size(0) != height) || (image_g.size(1) != width) ||
       (image_b.size(0) != height) || (image_b.size(1) != width))
      throw ex_invalid_argument("RGB image matrices dimensions mismatch");
   /* check arguments - quality */
   if ((quality < 0) || (quality > 1))
      throw ex_invalid_argument("image quality must be in [0,1]");
   /* open file for writing */
   FILE* f = fopen(filename, "w");
   if (f == NULL)
      throw ex_file_open_error(filename, "write");
   /* initialize jpeg compression */
   struct jpeg_compress_struct cinfo;
   struct jpeg_error_mgr jerr;
   cinfo.err = jpeg_std_error(&jerr);
   jpeg_init_compress(&cinfo);
   jpeg_stdio_dest(&cinfo, f);
   /* set image parameters */
   cinfo.image_height = height;
   cinfo.image_width  = width;
   cinfo.input_components = 3;
   cinfo.in_color_space = JCS_RGB;
   jpeg_set_defaults(&cinfo);
   jpeg_set_quality(&cinfo, static_cast<int>(quality*100), false);
   /* initialize line buffer */
   unsigned long buffer_size = 3*width; /* 3 components (RGB) per pixel */
   unsigned char* buffer = new unsigned char[buffer_size];
   /* write image data */
   jpeg_start_compress(&cinfo, true);
   for (unsigned long n = 0, x = 0; x < height; x++) {
      /* assemble scanline */
      for (unsigned long y = 0; y < buffer_size; n++) {
         double val_r = math::round(image_r[n] * UCHAR_MAX);
         double val_g = math::round(image_g[n] * UCHAR_MAX);
         double val_b = math::round(image_b[n] * UCHAR_MAX);
         buffer[y++] = static_cast<unsigned char>(val_r);
         buffer[y++] = static_cast<unsigned char>(val_g);
         buffer[y++] = static_cast<unsigned char>(val_b);
      }
      /* write scanline */
      if (jpeg_write_scanlines(&cinfo, &buffer, 1) != 1) {
         /* error writing image data */
         jpeg_finish_compress(&cinfo);
         jpeg_destroy_compress(&cinfo);
         delete [] buffer;
         fclose(f);
         throw ex_file_write_error(filename, "error writing jpeg image data");
      }
   }
   /* cleanup */
   jpeg_finish_compress(&cinfo);
   jpeg_destroy_compress(&cinfo);
   delete [] buffer;
   fclose(f);
}

} /* namespace image */
} /* namespace formats */
} /* namespace io */
