/*
 * Jpeg test.
 */
#include "io/formats/image/jpeg.hh"
#include "io/streams/cout.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/string.hh"
#include "math/matrices/matrix.hh"

using io::formats::image::jpeg;
using io::streams::cout;
using lang::exceptions::exception;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;
using lang::string;
using math::matrices::matrix;

int main(const int argc, const char** argv) {
   try {
      /* get filenames */
      if (argc < 3)
         throw ex_invalid_argument("must specify input and output filenames");
      string<> filename(argv[1]);
      string<> filename_out(argv[2]);
      string<> tab("   ");
      /* detect color/grayscale */
      cout << "processing: " << filename << "\n";
      bool is_gray = jpeg::is_grayscale(filename);
      cout << tab << ((is_gray) ? "grayscale image" : "color image") << "\n";
      /* process image */
      if (is_gray) {
         /* read grayscale image */
         auto_ptr< matrix<> > im_gray;
         jpeg::read(filename, im_gray);
         cout << tab << "size: "
              << im_gray->size(0) << " x " << im_gray->size(1) << "\n";
         /* write grayscale image */
         jpeg::write(filename_out, *im_gray);
         cout << tab << "wrote copy to: " << filename_out << "\n";
      } else {
         /* read color image */
         auto_ptr< matrix<> > im_r;
         auto_ptr< matrix<> > im_g;
         auto_ptr< matrix<> > im_b;
         jpeg::read(filename, im_r, im_g, im_b);
         cout << tab << "size: "
              << im_r->size(0) << " x " << im_r->size(1) << "\n";
         /* write color image */
         jpeg::write(filename_out, *im_r, *im_g, *im_b);
         cout << tab << "wrote copy to: " << filename_out << "\n";
      }
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
