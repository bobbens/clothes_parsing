/*
 * Serial input stream built on top of a standard istream.
 */
#include "io/serialization/serial_istream.hh"
#include "io/serialization/serial_ostream.hh"
#include "io/streams/istream.hh"
#include "lang/exceptions/ex_file_read_error.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::istream;
using lang::exceptions::ex_file_read_error;

namespace {
/*
 * Reverse endianness of the given data.
 */
template <typename T>
T reverse_endianness(const T& t) {
   T t_reverse;
   const char* t_ptr   = reinterpret_cast<const char*>(&t);
   char* t_reverse_ptr = reinterpret_cast<char*>(&t_reverse);
   for (unsigned long start = 0, end = sizeof(T); start < sizeof(T); start++) {
      t_reverse_ptr[start] = t_ptr[--end];
   }
   return t_reverse;
}

/*
 * Read data from the stream and reverse its endianness if needed.
 */
template <typename T>
void read_serial_istream(istream& is, T& t, bool reverse_endian) {
   is.read(reinterpret_cast<char*>(&t), sizeof(T));
   if (reverse_endian)
      t = reverse_endianness(t);
}
} /* namespace */

/*
 * Constructor.
 * Attach the istream for use as input.
 * The istream must be ready for reading.
 */
serial_istream::serial_istream(istream& is)
 : _is(is),
   _reverse_endian(false)
{
   /* read endian signature */
   unsigned short sig;
   _is.read(reinterpret_cast<char*>(&sig), sizeof(unsigned short));
   /* set endian reversal flag */
   unsigned short sig_correct = serial_istream::endian_signature();
   if (sig != sig_correct) {
      if (sig == reverse_endianness(sig_correct))
         _reverse_endian = true;
      else
         throw ex_file_read_error("invalid serial input stream");
   }
}

/*
 * Copy constructor.
 */
serial_istream::serial_istream(serial_istream& s)
 : _is(s._is),
   _reverse_endian(s._reverse_endian)
{ }

/*
 * Destructor.
 */
serial_istream::~serial_istream() {
   /* do nothing */
}

/*
 * Deserialization operators.
 * Read the specified data from the file stream.
 */
serial_istream& serial_istream::operator>>(bool& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(char& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(unsigned char& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(short& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(unsigned short& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(int& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(unsigned int& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(long& x) {
   long long y;
   *this >> y;
   x = static_cast<long>(y);
   return *this;
}
   
serial_istream& serial_istream::operator>>(unsigned long& x) {
   unsigned long long y;
   *this >> y;
   x = static_cast<unsigned long>(y);
   return *this;
}
   
serial_istream& serial_istream::operator>>(long long& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(unsigned long long& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(float& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(double& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

serial_istream& serial_istream::operator>>(long double& x) {
   read_serial_istream(_is, x, _reverse_endian);
   return *this;
}

/*
 * Return a value that changes under byte-reversal.
 * This returns the same value as serial_ostream::endian_signature().
 */
unsigned short serial_istream::endian_signature() {
   return serial_ostream::endian_signature();
}

} /* namespace serialization */
} /* namespace io */
