// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <assert.h>
extern "C"
{
#include <jpeglib.h>
#include <jerror.h>
}

#include "util.hh"
#include "image.hh"
#include "sort.hh"


namespace Util
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // read in a jpeg file into an ImageStack. 
  // if the jpeg file is grayscale then the resulting ImageStack only has 1 layer
  // otherwise it has 3 layers corresponding to the RGB colorspace.
  //

  bool readJpegFile(const char *filename, ImageStack& im)
  {
      assert(filename != NULL);

      int bytesPerPixel;
      int height;
      int width;
      unsigned char* imbuf;
      bool isGray;

      FILE* file = fopen(filename, "r");
      if (!file)
      {
        return (false);
      }

      jpeg_decompress_struct cinfo;
      jpeg_error_mgr jerr;
      cinfo.err = jpeg_std_error (&jerr);
      jpeg_create_decompress (&cinfo);
      jpeg_stdio_src (&cinfo, file);
      jpeg_read_header (&cinfo, (boolean) true);

      if (cinfo.out_color_space == JCS_GRAYSCALE)
      {
          isGray = true;
          cinfo.output_components = 1;
          bytesPerPixel = 1;
      }
      else
      {
          isGray = false;
          cinfo.out_color_space = JCS_RGB;
          cinfo.output_components = 3;
          bytesPerPixel = 3;
      }
      jpeg_calc_output_dimensions (&cinfo);
      jpeg_start_decompress (&cinfo);

      height = cinfo.output_height;
      width = cinfo.output_width;
      imbuf = new unsigned char[width * height * bytesPerPixel];

      const int lineSize = width * bytesPerPixel;
      unsigned char* p = imbuf;

      while (cinfo.output_scanline < cinfo.output_height)
      {
        jpeg_read_scanlines (&cinfo, &p, 1);
        p += lineSize;
      }

      jpeg_finish_decompress (&cinfo);
      jpeg_destroy_decompress (&cinfo);
      fclose (file);

      //now fill in our image array data structure
      if (isGray)
      {
        im.resize(1,width,height);
        for (int i = 0; i < width * height; i++)
        {
          const float gray = (float) imbuf[i] / 255;
          const int y = i / width;
          const int x = i % width;
          im(0,x,y) = gray;
          // this fails on fp type errors:  assert (gray >= 0 && gray <= 1);
        }
      }
      else
      {
        im.resize(3,width,height);
        for (int i = 0; i < width * height; i++)
        {
          const int y = i / width;
          const int x = i % width;
          float r = (float) imbuf[3 * i] / 255;
          float g = (float) imbuf[3 * i + 1] / 255;
          float b = (float) imbuf[3 * i + 2] / 255;
          if ((r < 0) || (r > 1))
          {
            std::cerr << "r = " << r << std::endl;
          }
          if ((g < 0) || (g > 1))
          {
            std::cerr << "g = " << g << std::endl;
          }
          if ((b < 0) || (b > 1))
          {
            std::cerr << "b = " << r << std::endl;
          }
          assert (r >= 0 && r <= 1);
          assert (g >= 0 && g <= 1);
          assert (b >= 0 && b <= 1);
          im(RGB_R,x,y) = r;
          im(RGB_G,x,y) = g;
          im(RGB_B,x,y) = b;
        }
      }
      delete[] imbuf;

      return true;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  //
  // used for writing out a single channel image with
  // the "jet" psuedo-color map
  //
  static int jetR (float v)
  {
    assert (v >= 0 && v <= 1);
    int i = (uint) rint (v * 255);
    return Util::minmax (0, (450 - 5 * abs (i - 196)), 255);
  }

  static int jetG (float v)
  {
    assert (v >= 0 && v <= 1);
    int i = (uint) rint (v * 255);
    return Util::minmax (0, (450 - 5 * abs (i - 128)), 255);
  }

  static int jetB (float v)
  {
      assert (v >= 0 && v <= 1);
      int i = (uint) rint (v * 255);
      return Util::minmax (0, (450 - 5 * abs (i - 64)), 255);
  }

  //
  // create a normalized version of this image.
  // first run thru and find the min and max values
  // and then create a new image whose values range 
  // over [0.0,1.0].  the constant image is set to 0.0
  //
  // used in the Jpeg reading and writing.
  //
  static void getNormalized(const Image& im, Image& normalized) 
  {
      int width = im.size(1);
      int height = im.size(0);
      normalized.resize(width,height);

      float max = im(0,0);
      float min = im(0,0);

      for (int i = 0; i < width; i++)
      {
        for (int j = 0; j < height; j++)
        {
          float val = im(i,j);
          max = Util::max(max, val);
          min = Util::min(min, val);
        }
      }

      if ((max - min) > 0)
      {
        normalized = (im - min) / (max - min);
      }
      else
      {
        normalized.init(0);
      }
  }

  //
  // write out a grayscale image to a jpeg file.
  // if normalize=true, then the range of the image is adjusted to use the full scale
  // if jet=true then the image is written in pseudocolor rather than grayscale
  //
  bool writeJpegFile(const Image& im, const char *filespec, const bool normalize, const bool jet)
  {
      assert(filespec != NULL);

      //normalized version of this image
      Image normim;
      if (normalize)
      {
        getNormalized(im,normim); 
      }
      else
      {
        normim = im;
      }

      struct jpeg_error_mgr jerr;
      struct jpeg_compress_struct cinfo;

      // Open the output file.
      FILE* file = Util::openOutputStrm (filespec);
      if (!file)
      {
        return (false);
      }

      // Set up the normal JPEG error handling.
      cinfo.err = jpeg_std_error (&jerr);

      // Init a JPEG compression object.
      jpeg_create_compress (&cinfo);

      // Specify source of data.
      jpeg_stdio_dest (&cinfo, file);

      // Specify compression parameters.
      cinfo.image_width = im.size(0);
      cinfo.image_height = im.size(1);
      cinfo.input_components = 3;
      cinfo.in_color_space = JCS_RGB;
      jpeg_set_defaults (&cinfo);

      // Start compression.
      jpeg_start_compress (&cinfo, TRUE);

      // Allocate space for one scanline.
      JSAMPLE* buf = new JSAMPLE[im.size(0) * 3];

      // Write data one scanline at a time.
      for (int y = 0; y < im.size(1); y++)
      {
        for (int x = 0; x < im.size(0); x++)
        {
          float v = normim(x,y);
          assert (v >= 0 && v <= 1);
          if (jet)
          {
              buf[3*x + 0] = jetR(v);
              buf[3*x + 1] = jetG(v);
              buf[3*x + 2] = jetB(v);
          }
          else
          {
              int g = (int)rint(255 * v);
              buf[3*x + 0] = g;
              buf[3*x + 1] = g;
              buf[3*x + 2] = g;
          }
        }
        int c = jpeg_write_scanlines (&cinfo, &buf, 1);
        if (c != 1)
        {
          fclose (file);
          jpeg_destroy_compress (&cinfo);
          delete[]buf;
          return false;
        }
      }

      // Clean up.
      jpeg_finish_compress (&cinfo);
      fclose (file);
      jpeg_destroy_compress (&cinfo);
      delete[]buf;
      return true;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // convert an RGB imagestack into an 1976 CIE L*a*b* imagestack.
  //

  void rgb2lab (const ImageStack& rgb, ImageStack& lab)
  {
      assert(rgb.size(0) == 3);
      const int w = rgb.size(1);
      const int h = rgb.size(2);

      lab.resize(3,w,h);

      for (int i = 0; i < w; i++)
      {
        for (int j = 0; j < h; j++)
        {
            // RGB
            const float R = rgb(RGB_R,i,j);
            const float G = rgb(RGB_G,i,j);
            const float B = rgb(RGB_B,i,j);
            assert (R >= 0 && R <= 1);
            assert (G >= 0 && G <= 1);
            assert (B >= 0 && B <= 1);
            // RGB --> XYZ
            const float X = 0.412453 * R + 0.357580 * G + 0.180423 * B;
            const float Y = 0.212671 * R + 0.715160 * G + 0.072169 * B;
            const float Z = 0.019334 * R + 0.119193 * G + 0.950227 * B;
            // XYZ of D65 reference white (R=G=B=1).
            const float Xn = 0.950456;
            const float Yn = 1.000000;
            const float Zn = 1.088754;
            // XYZ --> 1976 CIE L*a*b*
            const float rX = X / Xn;
            const float rY = Y / Yn;
            const float rZ = Z / Zn;
            const float thresh = 0.008856;
#define f(t) (t > thresh) ? pow(t,1./3.) : (7.787 * t + 16./116.)
            const float fX = f(rX);
            const float fY = f(rY);
            const float fZ = f(rZ);
#undef f
            const float L = (rY > thresh) ? 116. * pow (rY, 1. / 3.) - 16. : 903.3 * rY;
            const float a = 500. * (fX - fY);
            const float b = 200. * (fY - fZ);
            assert (L >= 0 && L <= 100);
            assert (a >= -120 && a <= 120);
            assert (b >= -120 && b <= 120);
            lab(LAB_L,i,j) = L;
            lab(LAB_A,i,j) = a;
            lab(LAB_B,i,j) = b;
        }
      }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // normalize an lab image so that values lie in [0,1]
  //
  void labNormalize (ImageStack& lab)
  {
    for(int x = 0; x < lab.size(1); x++)
    {
      for(int y = 0; y < lab.size(2); y++)
      {
        float L = lab(LAB_L,x,y);
        float a = lab(LAB_A,x,y);
        float b = lab(LAB_B,x,y);
        const float minab = -73;
        const float maxab = 95;
        const float range = maxab - minab;

        L = L / 100.0;
        a = (a - minab) / range;
        b = (b - minab) / range;
        L = Util::minmax(0.0f, L, 1.0f);
        a = Util::minmax(0.0f, a, 1.0f);
        b = Util::minmax(0.0f, b, 1.0f);
        lab(LAB_L,x,y) = L;
        lab(LAB_A,x,y) = a;
        lab(LAB_B,x,y) = b;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // create a translated version of this image where 
  // the old image appears embedded in a new image of
  // size [newwidth x newheight].  undefined pixels 
  // are filled in with value fill.
  //
  void getTranslated (const Image& im, const int xoffset, const int yoffset,
                      const int newwidth, const int newheight, const float fill, 
                      Image& translated)
  {

      assert(newwidth >= 0);
      assert(newheight >= 0);
      
      int oldwidth = im.size(1);
      int oldheight = im.size(0);
      translated.resize(newwidth,newheight);
      translated.init(fill);

      for (int x = xoffset; x < oldwidth + xoffset; x++)
      {
        for (int y = yoffset; y < oldheight + yoffset; y++)
        {
          if ((y >= 0) && (y < newheight) && (x >= 0) && (x < newwidth))
          {
            translated(x,y) = im(x-xoffset,y-yoffset);
          }
        }
      }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // utility method to compute an antialiasing filter of 
  // size 13 for the given downscaling factor
  //
  // used by getScaled
  //
  void createAntialiasFilter (const float scale, Image& filter) 
  {
      Util::Array1D<float> wind(13);
      Util::Array1D<float> b(6);
      float sum = 0;

      filter.resize(1,13);
      filter.init(0);

      for (int i = 0; i < 13; i++)
      {
          wind(i) = 0.54 - 0.46 * cos (2 * M_PI * i / 11);
      }

      for (int i = 0; i < 6; i++)
      {
          b(i) = sin (scale * M_PI * (i + 0.5)) / (M_PI * (i + 0.5));
      }

      for (int i = 0; i < 6; i++)
      {
          filter(1,i) = b(6 - i) * wind(i);
          filter(1,i + 6) = b(i) * wind(i + 6);
          sum = sum + filter(1,i) + filter(1,i + 6);
      }
      sum = fabs (sum);

      for (int i = 0; i < 12; i++)
      {
          filter(1,i) = filter(1,i) / sum;
      }
  }

  //
  // create a resized version of this image of size [newwidth x newheight]
  // if bilinear = true, use bilinear interpolation
  // otherwise use bicubic interpolation
  //
  void getScaled(const Image& im, const int newwidth, const int newheight, const bool bilinear, Image& scaled)
  {
      assert(newwidth >= 0);
      assert(newheight >= 0);

      //first compute the scale factor
      float oldheight = im.size(0);
      float oldwidth = im.size(1);
      float xscale = (float) newwidth / (float) oldwidth;
      float yscale = (float) newheight / (float) oldheight;

      //filter to prevent aliasing if necessary
      Image antialiased = im;
      
      if (xscale < 1)
      {
          Image filtx;
          Image filtered;
          createAntialiasFilter(xscale,filtx);
          getFiltered(antialiased,filtered,filtx);
          antialiased = filtered;
      }

      if (yscale < 1)
      {
          Image filty;
          Image filtered;
          createAntialiasFilter(yscale,filty);
          getFiltered(antialiased,filtered,filty);
          antialiased = filtered;
      }

      //build the affine matrix
      Matrix A(3,3);
      A.init(0);
      A(0,0) = xscale;
      A(1,1) = yscale;
      A(2,2) = 1.0;

      //transform the image
      getTransformed(antialiased,A,newwidth,newheight,0,0,bilinear,scaled);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // create a rotated version of this image
  // using an appropriate affine transform.
  // if bilinear = true, use bilinear interpolation
  // otherwise use bicubic.
  //
  void getRotated(const Image& im, const float angle, const bool bilinear, Image& rotated) 
  {
      //put theta in 0 - 2pi
      float theta = angle;
      while (theta < 0)
          theta = (theta + 2 * M_PI);
      while (theta > 2 * M_PI)
          theta = (theta - 2 * M_PI);

      //build the affine matrix
      Matrix A(3,3);
      A(0,0) =  cos(theta); A(0,1) = sin(theta); A(0,2) = 0.0;
      A(1,0) = -sin(theta); A(1,1) = cos(theta); A(1,2) = 0.0;
      A(2,0) = 0.0;         A(2,1) = 0.0;        A(2,2) = 1.0;

      getTransformed(im, A, im.size(0), im.size(1), 0, 0, bilinear, rotated);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // the cubic spline interpolation kernel
  //
  static float cubic_bspline (const float x)
  {
      float a, b, c, d;

      if ((x + 2.0) <= 0.0)
      {
          a = 0.0;
      }
      else
      {
          a = pow ((x + 2.0), 3.0);
      }

      if ((x + 1.0) <= 0.0)
      {
          b = 0.0;
      }
      else
      {
          b = pow ((x + 1.0), 3.0);
      }

      if (x <= 0)
      {
          c = 0.0;
      }
      else
      {
          c = pow (x, 3.0);
      }

      if ((x - 1.0) <= 0.0)
      {
          d = 0.0;
      }
      else
      {
          d = pow ((x - 1.0), 3.0);
      }

      return ((1.0 / 6.0) * (a - (4.0 * b) + (6.0 * c) - (4.0 * d)));
  }

  //
  // returns the inverse.  only works for 3x3 matricies
  //
  static void getInverse(const Matrix mat, Matrix& inv)
  {
      assert((mat.size(0) != 3) || (mat.size(1) != 3));
       
      inv.resize(3,3);

      float denom = mat(0, 0) * mat(1, 1) * mat(2, 2) -
                     mat(0, 0) * mat(1, 2) * mat(2, 1) -
                     mat(1, 0) * mat(0, 1) * mat(2, 2) +
                     mat(1, 0) * mat(0, 2) * mat(2, 1) +
                     mat(2, 0) * mat(0, 1) * mat(1, 2) -
                     mat(2, 0) * mat(0, 2) * mat(1, 1);
                                                                                                                         
      inv(0,0) = ( mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1)) / denom;
      inv(0,1) = (-mat(0, 1) * mat(2, 2) + mat(0, 2) * mat(2, 1)) / denom;
      inv(0,2) = ( mat(0, 1) * mat(1, 2) - mat(0, 2) * mat(1, 1)) / denom;
      inv(1,0) = (-mat(1, 0) * mat(2, 2) + mat(1, 2) * mat(2, 0)) / denom;
      inv(1,1) = ( mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0)) / denom;
      inv(1,2) = (-mat(0, 0) * mat(1, 2) + mat(0, 2) * mat(1, 0)) / denom;
      inv(2,0) = ( mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0)) / denom;
      inv(2,1) = (-mat(0, 0) * mat(2, 1) + mat(0, 1) * mat(2, 0)) / denom;
      inv(2,2) = ( mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0)) / denom;
  }

  //
  // returns a new image which is an affine
  // transformed version of this image.
  // newimage = A*image.  the new image
  // is of size (height, width) such that
  // the corners of the old image are transformed
  // to locations inside the new image
  //
  // if bilinear is TRUE then use bilinear interpolation
  // otherwise use bicubic B-spline interpolation
  //
  void getTransformed (const Image& im, const Matrix& A, const int width,
                         const int height, const int xoffset,
                         const int yoffset, const bool bilinear, 
                         Matrix& transformed) 
  {
      assert(width >= 0);
      assert(height >= 0);

      float oldwidth = (float)im.size(1);
      float oldheight = (float)im.size(0);

      Matrix B;  
      getInverse(A,B);

      //allocate the result
      transformed.resize(width,height);

      //transform the image
      for (float x = 0; x < width; x++)
      {
        for (float y = 0; y < height; y++)
        {
          //compute the coordinates in the original image plane
          float u = (x + xoffset) * B(0, 0) + (y + yoffset) * B(0, 1) + B(0, 2);
          float v = (x + xoffset) * B(1, 0) + (y + yoffset) * B(1, 1) + B(1, 2);

          //if it's outside the bounds of the 
          //source image, fill with zeros
          if ((u >= oldwidth) || (u < 0.0) || (v >= oldheight) || (v < 0.0))
          {
              transformed((int)x, (int)y) = 0.0;
          }
          else
          {
            //do bilinear or bicubic interpolation as required
            if (bilinear == true)
            {
              float u1 = floor (u);
              float u2 = u1 + 1;
              float v1 = floor (v);
              float v2 = v1 + 1;

              float du = u - u1;
              float dv = v - v1;

              u1 = Util::max (0.0f, u1);
              u2 = Util::min (oldwidth - 1, u2);
              v1 = Util::max (0.0f, v1);
              v2 = Util::min (oldheight - 1, v2);

              float val11 = im((int)u1, (int)v1);
              float val12 = im((int)u1, (int)v2);
              float val21 = im((int)u2, (int)v1);
              float val22 = im((int)u2, (int)v2);

              float val = (1 - dv) * (1 - du) * val11 + 
                            (1 - dv) * du * val12 + 
                             dv * (1 - du) * val21 + 
                              dv * du * val22;

              transformed((int)x, (int)y) = val;
            }
            else
            {
              float a = u - floor (u);
              float b = v - floor (v);
              float val = 0.0;
              for (int m = -1; m < 3; m++)
              {
                float r1 = cubic_bspline ((float) m - a);
                for (int n = -1; n < 3; n++)
                {
                  float r2 = cubic_bspline (-1.0 * ((float) n - b));
                  float u1 = floor (u) + m;
                  float v1 = floor (v) + n;
                  u1 = Util::min (oldwidth - 1, Util::max (0.0f, u1));
                  v1 = Util::min (oldheight - 1, Util::max (0.0f, v1));

                  val += im((int)u1,(int)v1) * r1 * r2;
                }
              }
              transformed((int)x,(int)y) = val;
            }
          }
        }
      }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // filters the image via convolution with the given 
  // kernel and returns the resulting image.  kernel 
  // must have odd dimensions.
  //
  void getFiltered (const Image& im, const Image& kernel, Image& filtered) 
  {
      // image and kernel dimensions
      const int iwidth = im.size(0);
      const int iheight = im.size(1);
      const int kwidth = kernel.size(0);
      const int kheight = kernel.size(1);

      //output is same size as input
      filtered.resize(iwidth,iheight);
      filtered.init(0);

      // the kernel must not be larger than the image
      assert (kwidth <= iwidth);
      assert (kheight <= iheight);

      // the kernel must be odd in each dimension so it can be centered
      // over a pixel
      assert ((kwidth % 2) == 1);
      assert ((kheight % 2) == 1);

      // radius of kernel in each dimension; also the coordinates of the
      // kernel's center pixel
      const int xr = kwidth / 2;
      const int yr = kheight / 2;

      //flip the kernel left-right and up-down
      Util::Array2D<float> kern(kwidth,kheight);
      kern.init(0);
      for (int x = 0; x < kwidth; x++)
      {
        const int xx = kwidth - 1 - x;
        for (int y = 0; y < kheight; y++)
        {
          const int yy = kheight - 1 - y;
          kern(xx,yy) = kernel(x,y);
        }
      }

      // padded image dimensions
      const int pwidth = iwidth + 2 * xr;
      const int pheight = iheight + 2 * yr;

      //create image with reflective padding
      Util::Array2D<float> pim(pwidth,pheight);

      // top left
      for (int y = 0; y < yr; y++)
      {
          const int py = yr - 1 - y;
          for (int x = 0; x < xr; x++)
            {
                const int px = xr - 1 - x;
                pim(px,py) = im(x,y);
            }
      }
      // top right
      for (int y = 0; y < yr; y++)
      {
          const int py = yr - 1 - y;
          for (int x = 0; x < xr; x++)
            {
                const int xs = iwidth - 1 - x;
                const int xd = xr + iwidth + x;
                pim(xd,py) = im(xs,y);
            }
      }
      // bottom left
      for (int y = 0; y < yr; y++)
      {
          const int ys = iheight - 1 - y;
          const int yd = yr + iheight + y;
          for (int x = 0; x < xr; x++)
            {
                const int px = xr - 1 - x;
                pim(px,yd) = im(x,ys);
            }
      }
      // bottom right
      for (int y = 0; y < yr; y++)
      {
          const int ys = iheight - 1 - y;
          const int yd = yr + iheight + y;
          for (int x = 0; x < xr; x++)
            {
                const int xs = iwidth - 1 - x;
                const int xd = xr + iwidth + x;
                pim(xd,yd) = im(xs,ys);
            }
      }
      // top
      for (int y = 0; y < yr; y++)
      {
          const int py = yr - 1 - y;
          for (int x = 0; x < iwidth; x++)
            {
                const int px = x + xr;
                pim(px,py) = im(x,y);
            }
      }
      // bottom
      for (int y = 0; y < yr; y++)
      {
          const int ys = iheight - 1 - y;
          const int yd = yr + iheight + y;
          for (int x = 0; x < iwidth; x++)
            {
                const int px = x + xr;
                pim(px,yd) = im(x,ys);
            }
      }
      // left 
      for (int y = 0; y < iheight; y++)
      {
          const int py = yr + y;
          for (int x = 0; x < xr; x++)
            {
                const int px = xr - 1 - x;
                pim(px,py) = im(x,y);
            }
      }
      // right
      for (int y = 0; y < iheight; y++)
      {
          const int py = yr + y;
          for (int x = 0; x < xr; x++)
            {
                const int xs = iwidth - 1 - x;
                const int xd = xr + iwidth + x;
                pim(xd,py) = im(xs,y);
            }
      }
      // center
      for (int y = 0; y < iheight; y++)
      {
          const int py = yr + y;
          for (int x = 0; x < iwidth; x++)
          {
              const int px = xr + x;
              pim(px,py) = im(x,y);
          }
      }


      // use direct access to underlying arrays for speed
      float *p_pim = pim.data();
      float *p_kern = kern.data();
      float *p_filtered = filtered.data();

      // do the convolution
      // interchange y and ky loops, and unroll ky loop
      // gets 371 MFLOPS (including all overhead) on 700MHz PIII (53% of peak)
      for (int x = 0; x < iwidth; x++)
      {
        for (int kx = 0; kx < kwidth; kx++)
        {
          const int pcol = (x + kx) * pheight;
          const int kcol = kx * kheight;
          int ky = 0;
          while (ky < (kheight & ~0x7))
          {
            const float k0 = p_kern[kcol + ky + 0];
            const float k1 = p_kern[kcol + ky + 1];
            const float k2 = p_kern[kcol + ky + 2];
            const float k3 = p_kern[kcol + ky + 3];
            const float k4 = p_kern[kcol + ky + 4];
            const float k5 = p_kern[kcol + ky + 5];
            const float k6 = p_kern[kcol + ky + 6];
            const float k7 = p_kern[kcol + ky + 7];
            float in0 = p_pim[pcol + 0 + ky + 0];
            float in1 = p_pim[pcol + 0 + ky + 1];
            float in2 = p_pim[pcol + 0 + ky + 2];
            float in3 = p_pim[pcol + 0 + ky + 3];
            float in4 = p_pim[pcol + 0 + ky + 4];
            float in5 = p_pim[pcol + 0 + ky + 5];
            float in6 = p_pim[pcol + 0 + ky + 6];
            for (int y = 0; y < iheight; y++)
            {
                const float in7 = p_pim[pcol + y + ky + 7];
                p_filtered[x*iheight + y] += k0*in0 + k1*in1 + k2*in2 + k3*in3 +
                                             k4*in4 + k5*in5 + k6*in6 + k7*in7;
                in0 = in1;
                in1 = in2;
                in2 = in3;
                in3 = in4;
                in4 = in5;
                in5 = in6;
                in6 = in7;
            }
            ky += 8;
          }
          while (ky < (kheight & ~0x3))
          {
            const float k0 = p_kern[kcol + ky + 0];
            const float k1 = p_kern[kcol + ky + 1];
            const float k2 = p_kern[kcol + ky + 2];
            const float k3 = p_kern[kcol + ky + 3];
            float in0 = p_pim[pcol + 0 + ky + 0];
            float in1 = p_pim[pcol + 0 + ky + 1];
            float in2 = p_pim[pcol + 0 + ky + 2];
            for (int y = 0; y < iheight; y++)
            {
                const float in3 = p_pim[pcol + y + ky + 3];
                p_filtered[x*iheight + y] += k0*in0 + k1*in1 + k2*in2 + k3*in3;
                in0 = in1;
                in1 = in2;
                in2 = in3;
            }
            ky += 4;
          }
          while (ky < (kheight & ~0x1))
          {
            const float k0 = p_kern[kcol + ky + 0];
            const float k1 = p_kern[kcol + ky + 1];
            float in0 = p_pim[pcol + 0 + ky + 0];
            for (int y = 0; y < iheight; y++)
            {
                const float in1 = p_pim[pcol + y + ky + 1];
                p_filtered[x*iheight + y] += k0*in0 + k1*in1;
                in0 = in1;
            }
            ky += 2;
          }
          while (ky < kheight)
          {
            const float k0 = p_kern[kcol + ky];
            for (int y = 0; y < iheight; y++)
            {
                p_filtered[x*iheight + y] += k0*p_pim[pcol + y + ky];
            }
            ky += 1;
          }
          assert (ky == kheight);
        }
      }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // filter at the given radius.
  // the resulting value at a pixel is the n'th-order value
  // in a circular window cetnered at the pixel
  // 0 <= order <= 1, 0 is smallest value, 1 is largest.
  // Call with 0.5 to get median filtering.
  //
  void getPctFiltered(const Image& im, const float radius, const float order, Image& filtered) 
  {
      assert (order >= 0 && order <= 1);

      if (radius < 1)
      {
        filtered = im;
        return;
      }

      //allocate window and filtered image
      const int windowRadius = (int) ceil (radius);
      const int windowDiam = 2 * windowRadius + 1;
      const int windowPixels = windowDiam * windowDiam;
      Util::Array1D<float> values(windowPixels);

      const int iwidth = im.size(0);
      const int iheight = im.size(1);
      filtered.resize(iwidth,iheight);

      //loop over the image
      for (int x = 0; x < iwidth; x++)
      {
        for (int y = 0; y < iheight; y++)
        {
          //copy values out of the window
          int count = 0;
          for (int u = -windowRadius; u <= windowRadius; u++)
          {
            for (int v = -windowRadius; v <= windowRadius; v++)
            {
              if ((u * u + v * v) > radius * radius)
              {
                continue;
              }
              int yi = y + u;
              int xi = x + v;
              if (yi < 0 || yi >= iheight)
              {
                continue;
              }
              if (xi < 0 || xi >= iwidth)
              {
                continue;
              }
              assert (count < windowPixels);
              values(count++) = im(xi,yi);
            }
          }
          //sort the values in ascending order
          assert (count > 0);
          Util::sort(values.data(), count);
          assert(values(0) <= values(count - 1));
          //pick out percentile value
          int index = (int) rint (order * (count - 1));
          assert (index >= 0 && index < count);
          float pctVal = values(index);
          assert (pctVal >= values(0));
          assert (pctVal <= values(count - 1));
          filtered(x,y) = pctVal;
        }
      }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // filter at the given radius.
  // the resulting value at a pixel is the maximum value
  // in a circular window cetnered at the pixel
  //
  void getMaxFiltered (const Image& im, const float radius, Image& filtered)
  {
    if (radius < 1)
    {
      filtered = im;
      return;
    }

    const int iwidth = im.size(0);
    const int iheight = im.size(1);
    const int windowRadius = (int)ceil(radius);
    filtered.resize(iwidth,iheight);

    //loop over the image
    for (int x = 0; x < iwidth; x++)
    {
      for (int y = 0; y < iheight; y++)
      {
        //extract max from window
        float maxVal = im(x,y);
        for (int u = -windowRadius; u <= windowRadius; u++)
        {
          for (int v = -windowRadius; v <= windowRadius; v++)
          {
            if ((u * u + v * v) > radius * radius)
            {
              continue;
            }
            int xi = x + u;
            int yi = y + v;
            if (yi < 0 || yi >= iheight)
            {
              continue;
            }
            if (xi < 0 || xi >= iwidth)
            {
              continue;
            }
            const float val = im(xi,yi);
            maxVal = Util::max(val,maxVal);
          }
        }
        // save max
        filtered(x,y) = maxVal;
      }
    }
  }

} //namespace Util

