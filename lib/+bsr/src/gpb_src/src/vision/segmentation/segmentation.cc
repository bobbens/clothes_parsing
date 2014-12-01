/*
 * Segmentation.
 */
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "mlearning/clustering/clusterers/graph/tree_clusterer.hh"
#include "vision/segmentation/boundary.hh"
#include "vision/segmentation/region.hh"
#include "vision/segmentation/segmentation.hh"

#include <iostream>

namespace vision {
namespace segmentation {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::array;
using lang::pointers::auto_ptr;
using math::matrices::matrix;
using mlearning::clustering::clusterers::graph::tree_clusterer;

/*
 * Compute region features given an oversegmentation.
 * Return the regions and boundaries.
 */
void compute_regions(
   const matrix<unsigned long>& assign,
   const matrix<>& im_L,
   const matrix<>& im_a,
   const matrix<>& im_b,
   const matrix<>& im_lc,
   const matrix<>& im_pb,
   auto_collection< region,   array_list<region> >&   regions,
   auto_collection< boundary, array_list<boundary> >& boundaries)
{
   /* initialize region and boundary arrays */
   regions.reset(new array_list<region>());
   boundaries.reset(new array_list<boundary>());
   /* allocate regions */
   unsigned long n_regions = max(assign) + 1;
   for (unsigned long n = 0; n < n_regions; n++) {
      auto_ptr<region> reg(new region);
      reg->_id = n;
      regions->add(*reg);
      reg.release();
   }
   /* get image size */
   unsigned long n_rows = assign.size(0);
   unsigned long n_cols = assign.size(1);
   /* loop over pixels, updating regions and boundaries */
   for (unsigned long r = 0; r < n_rows; r++) {
      for (unsigned long c = 0; c < n_cols; c++) {
         /* get pixel assignment */
         unsigned long pix_assign = assign(r,c);
         /* get the pixel image data */
         double pix_L = im_L(r,c);
         double pix_a = im_a(r,c);
         double pix_b = im_b(r,c);
         /* update associated region */
         region& reg = (*regions)[pix_assign];
         reg._size++;
         reg._sum_L  += pix_L;
         reg._sum_L2 += pix_L*pix_L;
         reg._sum_a  += pix_a;
         reg._sum_a2 += pix_a*pix_a;
         reg._sum_b  += pix_b;
         reg._sum_b2 += pix_b*pix_b;
         /* check boundary right */
         if (c < (n_cols-1)) {
            unsigned long pix_other_assign = assign(r,c+1);
            if (pix_other_assign != pix_assign) {
               /* get other region on boundary */
               region& rb = (*regions)[pix_other_assign];
               /* create boundary if needed */
               if (!reg._boundary_map.contains(pix_other_assign)) {
                  auto_ptr<boundary> b(new boundary);
                  b->_id = boundaries->size();
                  boundaries->add(*b);
                  reg._boundary_map.add(rb._id, b->_id);
                  rb._boundary_map.add(reg._id, b->_id);
                  b.release();
               }
               unsigned long& b_id = reg._boundary_map.find_image(pix_other_assign);
               boundary& b = (*boundaries)[b_id];
               /* update boundary */
               double pix_lc0 = im_lc(r,c);
               double pix_lc1 = im_lc(r,c+1);
               double pix_pb0 = im_pb(r,c);
               double pix_pb1 = im_pb(r,c+1);
               b._size++;
               b._sum_contrast += ((pix_lc0 > pix_lc1) ? pix_lc0 : pix_lc1);
               b._sum_pb       += ((pix_pb0 > pix_pb1) ? pix_pb0 : pix_pb1);
            }
         }
         /* check boundary down */
         if (r < (n_rows-1)) {
            unsigned long pix_other_assign = assign(r+1,c);
            if (pix_other_assign != pix_assign) {
               /* get other region on boundary */
               region& rb = (*regions)[pix_other_assign];
               /* create boundary if needed */
               if (!reg._boundary_map.contains(pix_other_assign)) {
                  auto_ptr<boundary> b(new boundary);
                  b->_id = boundaries->size();
                  boundaries->add(*b);
                  reg._boundary_map.add(rb._id, b->_id);
                  rb._boundary_map.add(reg._id, b->_id);
                  b.release();
               }
               unsigned long& b_id = reg._boundary_map.find_image(pix_other_assign);
               boundary& b = (*boundaries)[b_id];
               /* update boundary */
               double pix_lc0 = im_lc(r,c);
               double pix_lc1 = im_lc(r+1,c);
               double pix_pb0 = im_pb(r,c);
               double pix_pb1 = im_pb(r+1,c);
               b._size++;
               b._sum_contrast += ((pix_lc0 > pix_lc1) ? pix_lc0 : pix_lc1);
               b._sum_pb       += ((pix_pb0 > pix_pb1) ? pix_pb0 : pix_pb1);
            }
         }
         /* check boundary diagonal */
         if ((r < (n_rows-1)) && (c < (n_cols-1))) {
            unsigned long pix_other_assign = assign(r+1,c+1);
            if (pix_other_assign != pix_assign) {
               /* get other region on boundary */
               region& rb = (*regions)[pix_other_assign];
               /* create boundary if needed */
               if (!reg._boundary_map.contains(pix_other_assign)) {
                  auto_ptr<boundary> b(new boundary);
                  b->_id = boundaries->size();
                  boundaries->add(*b);
                  reg._boundary_map.add(rb._id, b->_id);
                  rb._boundary_map.add(reg._id, b->_id);
                  b.release();
               }
               unsigned long& b_id = reg._boundary_map.find_image(pix_other_assign);
               boundary& b = (*boundaries)[b_id];
               /* update boundary */
               double pix_lc0 = im_lc(r,c);
               double pix_lc1 = im_lc(r+1,c+1);
               double pix_pb0 = im_pb(r,c);
               double pix_pb1 = im_pb(r+1,c+1);
               b._size++;
               b._sum_contrast += ((pix_lc0 > pix_lc1) ? pix_lc0 : pix_lc1);
               b._sum_pb       += ((pix_pb0 > pix_pb1) ? pix_pb0 : pix_pb1);
            }
         }
      }
   }
}

/*
 * Region agglomerator.
 */
class region_agglm : public tree_clusterer<region,double>::agglomerator {
public:
   region_agglm(
      unsigned long n_regions,
      auto_collection< boundary, array_list<boundary> > boundaries)
    : _n_merges(0),
      _n_regions_init(n_regions),
      _n_regions(n_regions),
      _boundaries(boundaries)
   { }

   ~region_agglm() { }

   bool is_mergeable(const double& c) const {
      return ((_n_regions_init - _n_merges) > 10);
   }
      
   /*
    * Merge regions.
    */
   auto_ptr<region> merge(
      const tree_clusterer<region,double>::tree& t0,
      const tree_clusterer<region,double>::tree& t1) const
   {
      /* get regions to merge */
      region& r0 = *(t0.data);
      region& r1 = *(t1.data);
      /* merge regions */
      auto_ptr<region> r(new region);
      r->_id     = _n_regions++;
      r->_size   = r0._size + r1._size;
      r->_sum_L  = r0._sum_L  + r1._sum_L;
      r->_sum_L2 = r0._sum_L2 + r1._sum_L2;
      r->_sum_a  = r0._sum_a  + r1._sum_a;
      r->_sum_a2 = r0._sum_a2 + r1._sum_a2;
      r->_sum_b  = r0._sum_b  + r1._sum_b;
      r->_sum_b2 = r0._sum_b2 + r1._sum_b2;
      _n_merges++;
      return r;
   }

   /*
    * Compute boundary between regions.
    */
   auto_ptr<double> update(
      const tree_clusterer<region,double>::tree& t0,
      const tree_clusterer<region,double>::tree& t1) const
   {
      /* determine which is the newest vertex */
      const tree_clusterer<region,double>::tree& t_old = 
         (t0.data->_id < t1.data->_id) ? t0 : t1;
      const tree_clusterer<region,double>::tree& t_new =
         (t0.data->_id < t1.data->_id) ? t1 : t0;
      /* get regions to update boundary between */
      region& r_old = *(t_old.data);
      region& r_new = *(t_new.data);
      /* check if boundary needs to be computed */
      if (r_new._id >= _n_regions_init) {
         /* get subregions of new region */
         region& r_left  = *(t_new.left->data);
         region& r_right = *(t_new.right->data);
         /* check if subregions connected to old region */
         bool left_conn  = r_old._boundary_map.contains(r_left._id);
         bool right_conn = r_old._boundary_map.contains(r_right._id);
         if (left_conn && right_conn) {
            /* get the boundaries */
            unsigned long& b_left_id  = r_old._boundary_map.find_image(r_left._id);
            unsigned long& b_right_id = r_old._boundary_map.find_image(r_right._id);
            boundary& b_left  = (*_boundaries)[b_left_id];
            boundary& b_right = (*_boundaries)[b_right_id];
            /* merge the boundaries */
            auto_ptr<boundary> b(new boundary());
            b->_id = _boundaries->size();
            b->_size = b_left._size + b_right._size;
            b->_sum_contrast = b_left._sum_contrast + b_right._sum_contrast;
            b->_sum_pb = b_left._sum_pb + b_right._sum_pb;
            /* update connectivity */
            r_old._boundary_map.remove(r_left._id);
            r_old._boundary_map.remove(r_right._id);
            r_old._boundary_map.add(r_new._id, b->_id);
            r_new._boundary_map.add(r_old._id, b->_id);
            /* update boundary array */
            _boundaries->add(*b);
            b.release();
         } else if (left_conn) {
            /* boundary is left boundary */
            unsigned long& b_left_id = r_old._boundary_map.find_image(r_left._id);
            r_old._boundary_map.remove(r_left._id);
            r_old._boundary_map.add(r_new._id, b_left_id);
            r_new._boundary_map.add(r_old._id, b_left_id);
         } else if (right_conn) {
            /* boundary is right boundary */
            unsigned long& b_right_id = r_old._boundary_map.find_image(r_right._id);
            r_old._boundary_map.remove(r_right._id);
            r_old._boundary_map.add(r_new._id, b_right_id);
            r_new._boundary_map.add(r_old._id, b_right_id);
         }
      }
      /* lookup boundary */
      unsigned long& b_id = r_new._boundary_map.find_image(r_old._id);
      boundary& b = (*_boundaries)[b_id];
      /* compute cost */
      return region_agglm::cost(r_old, r_new, b);
   }
  
   /*
    * Compute region merge cost.
    */
   static auto_ptr<double> cost(
      const region& r0, const region& r1, const boundary& b)
   {
       
      static const double alpha_1 = 1;
      static const double alpha_2 = 1;
      static const double alpha_3 = 1;
      double mean_lc = b._sum_contrast/double(b._size);
      double mean_pb = b._sum_pb/double(b._size);
      double b_um = mean_lc + alpha_1 * mean_pb;
      double area0 = double(r0._size);
      double area1 = double(r1._size);
      double quad_err0
       = r0._sum_L2 - (r0._sum_L * r0._sum_L)/area0 + 
         r0._sum_a2 - (r0._sum_a * r0._sum_a)/area0 + 
         r0._sum_b2 - (r0._sum_b * r0._sum_b)/area0;
      double quad_err1
       = r1._sum_L2 - (r1._sum_L * r1._sum_L)/area1 + 
         r1._sum_a2 - (r1._sum_a * r1._sum_a)/area1 + 
         r1._sum_b2 - (r1._sum_b * r1._sum_b)/area1;
      double A0 = area0 + alpha_3 * quad_err0;
      double A1 = area1 + alpha_3 * quad_err1;
      double A = (A0 < A1) ? A0 : A1;
      double C = b_um * math::pow(A, alpha_2);
      return auto_ptr<double>(new double(C));
      
      /*
      static const double alpha_1 = 1;
      double mean_lc = b._sum_contrast/double(b._size);
      double mean_pb = b._sum_pb/double(b._size);
      double b_um = mean_lc + alpha_1 * mean_pb;
      return auto_ptr<double>(new double(b_um));
      */
      /*
      double n = r0._size + r1._size;
      double sum_L  = r0._sum_L  + r1._sum_L;
      double sum_L2 = r0._sum_L2 + r1._sum_L2;
      double sum_a  = r0._sum_a  + r1._sum_a;
      double sum_a2 = r0._sum_a2 + r1._sum_a2;
      double sum_b  = r0._sum_b  + r1._sum_b;
      double sum_b2 = r0._sum_b2 + r1._sum_b2;
      double C = sum_L2/n - (sum_L*sum_L)/(n*n)
               + sum_a2/n - (sum_a*sum_a)/(n*n)
               + sum_b2/n - (sum_b*sum_b)/(n*n);
      return auto_ptr<double>(new double(C));
      */
   }
protected:
   mutable unsigned long _n_merges;
   unsigned long _n_regions_init;
   mutable unsigned long _n_regions;
   mutable auto_collection< boundary, array_list<boundary> > _boundaries;
};

/*
 * Create a segmentation given initial regions and boundaries.
 * Return the assignment of initial regions -> final regions.
 */
array<unsigned long> segment(
   auto_collection< region,   array_list<region> >   regions,
   auto_collection< boundary, array_list<boundary> > boundaries)
{
   /* build edge connectivity */
std::cout << "init edges\n";
   auto_collection<
      array<unsigned long>, array_list< array<unsigned long> > > edges(
         new array_list< array<unsigned long> >()
      );
std::cout << "building edges\n";
   for (array_list<region>::iterator_t i(*regions); i.has_next(); ) {
      region& r = i.next();
      auto_ptr< array<unsigned long> > edge_inds(
         new array<unsigned long>(r._boundary_map.size())
      );
      map<unsigned long, unsigned long>::iterator_t ii(r._boundary_map);
      unsigned int n = 0;
      while (ii.has_next()) {
         (*edge_inds)[n] = ii.next();
         n++;
      }
      edges->add(*edge_inds);
      edge_inds.release();
   }
   /* cluster */
std::cout << "creating agglm\n";
   region_agglm r_agglm(regions->size(), boundaries);
std::cout << "calling tree clusterer\n";
   tree_clusterer<region,double> c(r_agglm);
   return c.cluster(*regions, *edges);
}

} /* namespace segmentation */
} /* namespace vision */
