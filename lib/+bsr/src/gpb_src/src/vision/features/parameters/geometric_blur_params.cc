/*
 * Geometric blur parameters.
 */
#include "lang/array.hh"
#include "vision/features/parameters/geometric_blur_params.hh"

namespace vision {
namespace features {
namespace parameters {
/*
 * Imports.
 */
using lang::array;

/*
 * Constructor.
 * Return the default set of parameters.
 */
geometric_blur_params::geometric_blur_params()
 : _alpha(0.5),         /* default blur amount at distance |x| to be */
   _beta(1.0),          /* 0.5*|x| + 1                               */
   _radii(),
   _oris(),
   _channels(4),        /* assume 4 signal channels */
   _descriptor_size(0)
{
   /* set radii (in pixels) at which geometric blur is sampled */
   /* and number of orientations sampled at each radius        */
   const unsigned long n_radii = 6;
   const double        r[n_radii] = { 0, 5.6, 11.2, 22.4, 44.8, 70 };
   const unsigned long o[n_radii] = { 1, 8,    8,   10,   12,   12 };
   _radii = array<double>(n_radii);
   _oris  = array<unsigned long>(n_radii);
   for (unsigned long n = 0; n < n_radii; n++) {
      _radii[n] = r[n];
      _oris[n]  = o[n];
   }
   /* compute number of bins in feature vector */
   this->compute_descriptor_size();
}

/*
 * Constructor.
 * Specify parameters.
 */
geometric_blur_params::geometric_blur_params(
   double                      alpha,
   double                      beta,
   const array<double>&        radii,
   const array<unsigned long>& oris,
   unsigned long               channels) 
 : _alpha(alpha), 
   _beta(beta), 
   _radii(radii), 
   _oris(oris),
   _channels(channels)
{
   /* compute number of bins in feature vector */
   this->compute_descriptor_size();
}

/*
 * Compute feature vector size from radii and orientation arrays.
 */
void geometric_blur_params::compute_descriptor_size() {
   unsigned long n_radii = _radii.size();
   _descriptor_size = 0;
   for (unsigned long n = 0; n < n_radii; n++)
      _descriptor_size += _oris[n];
   _descriptor_size *= _channels;
}

/*
 * Copy constructor.
 */
geometric_blur_params::geometric_blur_params(const geometric_blur_params& params) 
 : _alpha(params._alpha), 
   _beta(params._beta), 
   _radii(params._radii), 
   _oris(params._oris),
   _channels(params._channels),
   _descriptor_size(params._descriptor_size)
{ }

/*
 * Destructor.
 */
geometric_blur_params::~geometric_blur_params() {
   /* do nothing */
}

/*
 * Get parameters used for controlling variation of blur with distance.
 * The blur amount at distance |x| is alpha*|x| + beta.
 */
double geometric_blur_params::alpha() const {
   return _alpha;
}

double geometric_blur_params::beta() const {
   return _beta;
}

/*
 * Get the set of radii at which geometric blur is sampled.
 */
const array<double>& geometric_blur_params::radii() const {
   return _radii;
}

/*
 * Get the number of orientations to sample at each radius.
 */
const array<unsigned long>& geometric_blur_params::orientations() const {
   return _oris;
}

/*
 * Get number of signal channels.
 */
unsigned long geometric_blur_params::channels() const {
   return _channels;
}

/*
 * Return the number of bins in the geometric blur feature vector.
 */
unsigned long geometric_blur_params::descriptor_size() const {
   return _descriptor_size;
}

/*
 * Compute blur amount at the given distance from the feature center.
 */
double geometric_blur_params::blur_sigma(double x) const {
   return (_alpha*x + _beta);
}

/*
 * Default parameter set.
 */
const geometric_blur_params geometric_blur_params::default_parameters;

} /* namespace parameters */
} /* namespace features */
} /* namespace vision */
