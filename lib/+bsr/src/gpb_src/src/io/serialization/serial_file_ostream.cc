/*
 * Serial file output stream.
 */
#include "io/serialization/serial_file_ostream.hh"
#include "io/serialization/serial_ostream.hh"
#include "io/streams/ofstream.hh"
#include "lang/exceptions/ex_file_open_error.hh"
#include "lang/string.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ofstream;
using lang::exceptions::ex_file_open_error;
using lang::string;

/***************************************************************************
 * Implementation of serial_file_ostream_base.
 ***************************************************************************/
/*
 * Constructor.
 * Open the specified file for output.
 */
serial_file_ostream_base::serial_file_ostream_base(const string<>& filename)
 : _f(filename, ofstream::out | ofstream::binary)
{
   if (_f.fail())
      throw ex_file_open_error(filename, "write");
}

/*
 * Destructor.
 * Close the file.
 */
serial_file_ostream_base::~serial_file_ostream_base() {
   _f.close();
}

/***************************************************************************
 * Implementation of serial_file_ostream.
 ***************************************************************************/

/*
 * Constructor.
 * Open the specified file for use in serialization.
 */
serial_file_ostream::serial_file_ostream(const string<>& filename)
 : serial_file_ostream_base(filename),
   serial_ostream(this->_f)
{ }

/* 
 * Destructor.
 * Close the file.
 */
serial_file_ostream::~serial_file_ostream() {
   /* do nothing (file is closed by base class destructor) */
}

/*
 * Serialization operators.
 * Write the specified data to the file stream.
 */
serial_file_ostream& serial_file_ostream::operator<<(const bool& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const char& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const unsigned char& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const short& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const unsigned short& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const int& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const unsigned int& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const long& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const unsigned long& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const long long& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const unsigned long long& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const float& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const double& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

serial_file_ostream& serial_file_ostream::operator<<(const long double& x) {
   this->serial_ostream::operator<<(x);
   return *this;
}

} /* namespace serialization */
} /* namespace io */
