/*
 * Serial output stream built on top of a standard ostream.
 */
#include "io/serialization/serial_ostream.hh"
#include "io/streams/ostream.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Constructor.
 * Attach the ostream for use as output.
 * The ostream must be ready for writing.
 */
serial_ostream::serial_ostream(ostream& os)
 : _os(os)
{
   /* write endian signature */
   unsigned short sig = serial_ostream::endian_signature();
   _os.write(reinterpret_cast<const char*>(&sig), sizeof(unsigned short));
}

/*
 * Copy constructor.
 */
serial_ostream::serial_ostream(serial_ostream& s)
 : _os(s._os)
{ }

/*
 * Destructor.
 */
serial_ostream::~serial_ostream() {
   /* do nothing */
}

/*
 * Serialization operators.
 * Write the specified data to the stream.
 */
serial_ostream& serial_ostream::operator<<(const bool& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(bool));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const char& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(char));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const unsigned char& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(unsigned char));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const short& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(short));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const unsigned short& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(unsigned short));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const int& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(int));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const unsigned int& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(unsigned int));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const long& x) {
   return (*this << static_cast<long long>(x));
}

serial_ostream& serial_ostream::operator<<(const unsigned long& x) {
   return (*this << static_cast<unsigned long long>(x));
}

serial_ostream& serial_ostream::operator<<(const long long& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(long long));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const unsigned long long& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(unsigned long long));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const float& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(float));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const double& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(double));
   return *this;
}

serial_ostream& serial_ostream::operator<<(const long double& x) {
   _os.write(reinterpret_cast<const char*>(&x), sizeof(long double));
   return *this;
}

/*
 * Return a value that changes under byte-reversal.
 * This returns the same value as serial_istream::endian_signature().
 */
unsigned short serial_ostream::endian_signature() {
   unsigned short sig = 0xFF00;
   return sig;
}

} /* namespace serialization */
} /* namespace io */
