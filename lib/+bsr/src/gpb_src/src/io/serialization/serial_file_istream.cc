/*
 * Serial file input stream.
 */
#include "io/serialization/serial_file_istream.hh"
#include "io/serialization/serial_istream.hh"
#include "io/streams/ifstream.hh"
#include "lang/exceptions/ex_file_open_error.hh"
#include "lang/string.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ifstream;
using lang::exceptions::ex_file_open_error;
using lang::string;

/***************************************************************************
 * Implementation of serial_file_istream_base.
 ***************************************************************************/

/*
 * Constructor.
 * Open the specified file for input.
 */
serial_file_istream_base::serial_file_istream_base(const string<>& filename)
 : _f(filename, ifstream::in | ifstream::binary)
{
   if (_f.fail())
      throw ex_file_open_error(filename, "read");
}

/*
 * Destructor.
 * Close the file.
 */
serial_file_istream_base::~serial_file_istream_base() {
   _f.close();
}

/***************************************************************************
 * Implementation of serial_file_istream.
 ***************************************************************************/

/*
 * Constructor.
 * Open the specified file for use in deserialization.
 */
serial_file_istream::serial_file_istream(const string<>& filename)
 : serial_file_istream_base(filename),
   serial_istream(this->_f)
{ }

/* 
 * Destructor.
 * Close the file.
 */
serial_file_istream::~serial_file_istream() {
   /* do nothing (file is closed by base class destructor) */
}

/*
 * Deserialization operators.
 * Read the specified data from the file stream.
 */
serial_file_istream& serial_file_istream::operator>>(bool& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(char& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(unsigned char& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(short& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(unsigned short& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(int& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(unsigned int& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(long& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(unsigned long& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(long long& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(unsigned long long& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(float& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(double& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

serial_file_istream& serial_file_istream::operator>>(long double& x) {
   this->serial_istream::operator>>(x);
   return *this;
}

} /* namespace serialization */
} /* namespace io */
