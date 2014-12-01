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

#include <string.h>
#include <iostream>
#include "segmentation.hh"
#include "string.hh"
#include "exception.hh"
#include "util.hh"
#include "regex.hh"

namespace Util
{
  //
  // create an empty segmentation 
  // 
  Segmentation::Segmentation()
  {
    _map.resize(0,0);
    _map.init(0);

    _date = NULL;
    _image = NULL;
    _user = NULL;
    _imageID = 0; 
    _userID = 0;
    _presentation = 0;
  }

  //
  // create an empty segmentation of given size
  Segmentation::Segmentation(const int width, const int height)
  {
    _map.resize(width,height);
    _map.init(0);

    _date = NULL;
    _image = NULL;
    _user = NULL;
    _imageID = 0; 
    _userID = 0;
    _presentation = 0;

    _renumber();
  }

  //
  // create a segmentation which uses the
  // given partition as a map
  //
  Segmentation::Segmentation(const Array2D<int> partition)
  {
    _map = partition;

    _date = NULL;
    _image = NULL;
    _user = NULL;
    _imageID = 0; 
    _userID = 0;
    _presentation = 0;

    _renumber();
  }

  //
  // cleanup strings if necessary
  //
  Segmentation::~Segmentation ()
  {
    if (_date) { free (_date); }
    if (_image) { free (_image); }
    if (_user) { free (_user); }
  }


  // Put single-line string into a canonical form that is good for
  // pattern-matching.  Caller must free() returned string.
  static char* 
  _cleanupLine (const char* line)
  {
      static Regex re1 ("#.*$");
      static Regex re2 ("[ \t\r\f\v\n]+");
      static Regex re3 ("^ ");
      static Regex re4 (" $");
      char* s = strdup (line);
      re1.subst (s, "", false); 	// Remove comment.
      re2.subst (s, " ", true);	// Compress whitespace.
      re3.subst (s, "", false);	// Remove leading whitespace.
      re4.subst (s, "", false);	// Remove trailing whitespace.
      return s;
  }

  //
  // read in a segfile 
  //
  void Segmentation::readFile(const char* file)
  {
      String buf;
      int lineNo = 0;
      int width = 0;
      int height = 0;

      FILE* strm = Util::openInputStrm (file);

      // Parse file header.
      Regex reBlank 	("^$");
      Regex reDate 	("^date (.*)$");
      Regex reImage 	("^image (.*)$");
      Regex reUser 	("^user (.*)$");
      Regex reSegments 	("^segments ([0-9]+)$");
      Regex reWidth 	("^width ([0-9]+)$");
      Regex reHeight 	("^height ([0-9]+)$");
      Regex reGray 	("^gray (0|1)$");
      Regex reInvert 	("^invert (0|1)$");
      Regex reFlipFlop	("^flipflop (0|1)$");
      Regex reFormat 	("^format (ascii|binary) (cr|map)$");
      Regex reData 	("^data$");

      int format = defaultFormat;
      int k = 0;
      while (buf.nextLine (strm)) 
      {
        lineNo++;
        char* s = _cleanupLine (buf.text());
        char* sub = strdup (s);		// Scratch array for submatches.
        if (reBlank.match (s)) 
        {
            free (s); free (sub);
            continue;
        } 
        else if (reDate.match (s)) 
        {
            _date = strdup (reDate.get (1, sub));
        } 
        else if (reImage.match (s)) 
        {
            _image = strdup (reImage.get (1, sub));
            _imageID = atoi (_image);
        } 
        else if (reUser.match (s)) 
        {
            _user = strdup (reUser.get (1, sub));
            _userID = atoi (_user);
        } 
        else if (reSegments.match (s)) 
        {
            k = atoi (reSegments.get (1, sub));
        } 
        else if (reWidth.match (s)) 
        {
            width = atoi (reWidth.get (1, sub));
        } 
        else if (reHeight.match (s)) 
        {
            height = atoi (reHeight.get (1, sub));
        } 
        else if (reGray.match (s)) 
        {
            int x = atoi (reGray.get (1, sub));
            if (x) { _presentation |= presGray; } 
            else { _presentation &= ~presGray; }
        } 
        else if (reInvert.match (s)) 
        {
            int x = atoi (reInvert.get (1, sub));
            if (x) { _presentation |= presInvert; } 
            else { _presentation &= ~presInvert; }
        } 
        else if (reFlipFlop.match (s)) 
        {
            int x = atoi (reFlipFlop.get (1, sub));
            if (x) { _presentation |= presFlipFlop; } 
            else { _presentation &= ~presFlipFlop; }
        } 
        else if (reFormat.match (s)) 
        {
            format = 0;
            if (strcmp (reFormat.get(1,sub), "ascii") == 0) {
              format |= ffASCII;
            } else {
              format &= ~ffASCII;
            }
            if (strcmp (reFormat.get(2,sub), "cr") == 0) {
              format |= ffCr;
            } else {
              format &= ~ffCr;
            }
        } 
        else if (reData.match (s)) 
        {
            free (s);
            free (sub);
            break;
        } 
        else 
        {
            free (s);
            free (sub);
            throw Exception (String ("Segmentation::readFile: Error parsing header (line %d): %s.",lineNo, buf.text()));
        }
        free (s); 
        free (sub);
      }

      // Make sure we parsed out the critical info from the header.
      if (width == 0 || height == 0) 
      {
        throw Exception ("Segmentation::readFile: Image dimensions not specified."); 
      }
      if (k == 0) 
      {
        throw Exception ("Segmentation::readFile: Number of segments not specified.");
      }

      // We only support ascii/cr format.
      if (format != (ffASCII | ffCr)) {
        throw Exception ("Segmentation::readFile: Only supported .seg format is ascii/cr.");
      }

      //initialize segment map
      _map.resize(width,height);
      _map.init(k);

      // Read the data.
      Regex reDatum ("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)$");
      assert (format == (ffASCII | ffCr));
      while (buf.nextLine (strm)) 
      {
        lineNo++;
        char* s = _cleanupLine (buf.text());
        char* sub = strdup (s);		// Scratch array for submatches.
        if (strcmp (s, "") == 0) 
        {
          free (s); free (sub);
          continue;
        }
        if (!reDatum.match (s)) 
        {
          free (s); free (sub);
          throw Exception (String ("Segmentation::_init: Parse error (line %d): %s.",lineNo, buf.text()));
        }
        int seg = atoi (reDatum.get (1, sub));
        int y = atoi (reDatum.get (2, sub));
        int x1 = atoi (reDatum.get (3, sub));
        int x2 = atoi (reDatum.get (4, sub));
        free (s); free (sub);

        if (seg >= k) 
        {
          throw Exception (String (
          "Segmentation::_init: Line %d: %d > num segments (%d).", 
          lineNo, seg, k));
        }
        if (y >= height || x1 >= width || x2 >= width) 
        {
          throw Exception (String (
          "Segmentation::_init: Line %d: pixels [%d..%d,%d] "
          "is out of range [%d,%d].",
          lineNo, x1, x2, y, width, height));
        }

        for (int x = x1; x <= x2; x++) 
        {
          if (_map(x,y) < k) 
          {
            throw Exception (String (
              "Segmentation::_init: Line %d: pixel [%d,%d] "
              "assigned to 2 segments.",
              lineNo, x, y));
          }
          _map(x,y) = seg;
        }
      }    

      // Check that all pixels were assigned to a segment.
      for (int x = 0; x < width; x++) 
      {
        for (int y = 0; y < height; y++) 
        {
          if (_map(x,y) == k) 
          {
            throw Exception (String (
            "Segmentation::_init: pixel [%d,%d] not assigned "
            "to a segment.", x,y));
          }
        }
      }
      fclose (strm);
      _renumber();
  }


  //
  // write out to a segfile
  //
  void 
  Segmentation::writeFile (const char* file, int format) const
  {
    FILE* strm = Util::openOutputStrm (file);

    if (format != (ffASCII | ffCr)) // Check format.
    {
      throw Exception ("Segmentation::writeFile: Only supported .seg format is ascii/cr.");
    }
    assert (format == (ffASCII | ffCr));
    fprintf (strm, "format ascii cr\n");
    if (_date != NULL && strcmp (_date, "") != 0) 
    { 
      fprintf (strm, "date %s\n", _date); 
    }
    if (_image != NULL) { fprintf (strm, "image %s\n", _image); } 
    else if (_imageID > 0) { fprintf (strm, "image %d\n", _imageID); }
    if (_user != NULL) { fprintf (strm, "user %s\n", _user); }
    else if (_userID > 0) { fprintf (strm, "user %d\n", _userID); }
    fprintf (strm, "width %d\n", _map.size(0));
    fprintf (strm, "height %d\n", _map.size(1));
    fprintf (strm, "segments %d\n", _numSegments);
    fprintf (strm, "gray %d\n", isGray() ? 1 : 0);
    fprintf (strm, "invert %d\n", isInverted() ? 1 : 0);
    fprintf (strm, "flipflop %d\n", isFlipFloped() ? 1 : 0);
    fprintf (strm, "data\n");
    //write out in scanline order
    for (int row = 0; row < _map.size(1); row++) 
    {
      int prevSeg = _map(0,row);
      int start = 0;
      for (int col = 1; col < _map.size(0); col++) 
      {
        int seg = _map(col,row);
        if (seg == prevSeg) { continue; }
        fprintf (strm, "%u %u %u %u\n", prevSeg, row, start, col - 1);
        start = col;
        prevSeg = seg;
      }
      fprintf (strm, "%u %u %u %u\n", prevSeg, row, start, _map.size(0)-1);
    }
    fclose (strm);
  }


  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // connectedness and merging small segments
  //

  class Node {
    public:
      Node (int x, int y, Node* next) 
      {
        _x = x;
        _y = y;
        _next = next;
      }
      int _x;
      int _y;
      Node* _next;
  };

  //
  // breadth-first flood fill starting at (x,y), replacing oldId
  // with newId.
  //
  void
  Segmentation::_replace(int x, int y, int oldId, int newId, float radius)
  {
      assert (x < _map.size(0));
      assert (y < _map.size(1));
      assert (_map(x,y) == oldId);

      Node* oldList = new Node (x, y, NULL);
      _map(x,y) = newId;

      int window = (int) ceil (radius);
      float radiusSqr = radius * radius;
      while (oldList != NULL) 
      {
        Node* newList = NULL;
        while (oldList != NULL) 
        {
          //pop a node
          Node* node = oldList;
          int x = node->_x;
          int y = node->_y;
          oldList = node->_next;
          delete node;
          assert (x < _map.size(0));
          assert (y < _map.size(1));
          assert (_map(x,y) == newId);

          //push the node's neighbors
          int umin = Util::max (x - window, 0);
          int umax = Util::min (x + window, (int)_map.size(0) - 1);
          int vmin = Util::max (y - window, 0);
          int vmax = Util::min (y + window, (int)_map.size(1) - 1);
          for (int u = umin; u <= umax; u++) 
          {
            for (int v = vmin; v <= vmax; v++) 
            {
              float distSqr = (x-u)*(x-u) + (y-v)*(y-v);
              if (distSqr > radiusSqr) { continue; }
              if (_map(u,v) != oldId) { continue; }
              newList = new Node(u,v,newList);
              _map(u,v) = newId;
            }
          }
        }
        oldList = newList;
      }
  }

  //
  // renumber the segments in a contigous manner
  // remove size 0 segments, etc
  //
  void 
  Segmentation::_renumber()
  {
    //upper bound the number of segments
    int maxseg = 0;
    for (int x = 0; x < _map.size(0); x++)
    {
      for (int y = 0; y < _map.size(1); y++)
      {
        maxseg = Util::max(maxseg,_map(x,y));
      } 
    }

    //see which segment numbers are in use
    Array1D<bool> empty(maxseg+1);
    empty.init(true);
    for (int x = 0; x < _map.size(0); x++)
    {
      for (int y = 0; y < _map.size(1); y++)
      {
        empty(_map(x,y)) = false;
      } 
    }


    //compute renumbering from old segment numbers to new
    _numSegments = 0;
    Array1D<int> renum(maxseg+1);
    renum.init(-1);
    for (int s = 0; s < maxseg+1; s++)
    {
      if (!empty(s))
      {
        renum(s) = _numSegments;
        _numSegments++;
      }
    }

    //apply renumbering
    for (int x = 0; x < _map.size(0); x++)
    {
      for (int y = 0; y < _map.size(1); y++)
      {
        _map(x,y) = renum(_map(x,y));
        assert(_map(x,y) >= 0);
        assert(_map(x,y) < _numSegments);
      }
    }
  }


  //
  // break a segmentation into connected components so that
  // every segment label refers to exactly one connected 
  // component.
  //
  void
  Segmentation::fragment (float radius)
  {
    int nextId = _numSegments;
    for (int x = 0; x < _map.size(0); x++) 
    {
      for (int y = 0; y < _map.size(1); y++) 
      {
        int oldId = _map(x,y);
        if (oldId < _numSegments)
        {
          _replace(x,y,oldId,nextId,radius);
          nextId++;
        }
      }
    }
    _renumber();
  }

  //
  // merge segments containing fewer than dustthresh pixels
  // into neighboring regions
  //
  void
  Segmentation::mergeDust(int dustthresh)
  {
    fragment(1.0f);

    Array1D<int> sizes;
    segSizes(sizes);

    //if segment s is dustlike, merge with a neighbor
    Util::Message::debug(String("merging dust starting with %d segments",_numSegments),1);

    Array1D<int> nbrcnt(_numSegments);
    for (int s = 0; s < _numSegments; s++)
    {
      if (sizes(s) < dustthresh)
      {
        nbrcnt.init(0);

        //find how many nbring pixels this pixel has from each segment
        for (int x = 1; x < _map.size(0)-1; x++)
        {
          for (int y = 1; y < _map.size(1)-1; y++)
          {
            //if we are a neigbhor of s, add 1 to our nbrcnt
            if ( (_map(x  ,y+1) == s) || (_map(x+1,y+1) == s) ||
                 (_map(x+1,y  ) == s) || (_map(x+1,y-1) == s) ||
                 (_map(x  ,y-1) == s) || (_map(x-1,y-1) == s) ||
                 (_map(x-1,y  ) == s) || (_map(x-1,y+1) == s) )
            {
              nbrcnt(_map(x,y))++;   
            }
          }
        }

        //pick the most commonly occuring neighbor
        int bigestnbr = 0;
        for (int ss = 0; ss < _numSegments; ss++)
        {
          if (nbrcnt(ss) > nbrcnt(bigestnbr))
          {
            bigestnbr = ss;
          }
        }
        
        //and merge with that nbr
        for (int x = 0; x < _map.size(0); x++)
        {
          for (int y = 0; y < _map.size(1); y++)
          {
            if (_map(x,y) == s)
            {
              _map(x,y) = bigestnbr;
            }
          }
        }
      }
    }

    //lastly, renumber
    _renumber();
    Util::Message::debug(String("merged down to %d segments",_numSegments),1);
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  // segment -> boundary conversions
  //


  void Segmentation::computeBoundaryMap(Util::Array2D<bool>& boundaryMap) const
  {
    const int w = getWidth();
    const int h = getHeight();
    boundaryMap.resize(w,h);
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            const int p = _map(x,y);
            if ((x < w-1 && p != _map(x+1,y)) ||
                (y < h-1 && p != _map(x,y+1)) ||
                (x < w-1 && y < h-1 && p != _map(x+1,y+1)))
            {
                boundaryMap(x,y) = true;
            }
            else
            {
                boundaryMap(x,y) = false;
            }
        }
    }
  }

  void Segmentation::computeBoundaryMapHalf(Util::Array2D<bool>& boundaryMap) const
  {
    const int w = getWidth();
    const int h = getHeight();
    boundaryMap.resize(w/2,h/2);
    for (int x = 0; x/2 < w/2; x+=2) {
        for (int y = 0; y/2 < h/2; y+=2) {
            const int p = _map(x,y);
            if ((x < w-1 && p != _map(x+1,y)) ||
                (y < h-1 && p != _map(x,y+1)) ||
                (x < w-1 && y < h-1 && p != _map(x+1,y+1)))
            {
                boundaryMap(x/2,y/2) = true;
            }
            else
            {
                boundaryMap(x/2,y/2) = false;
            }
        }
    }
  }

  class Point {
    public:
      Point () { x = y = -1; }
      Point (int x, int y) { this->x = x; this->y = y; }
      int x, y;
  };

  class DistNode 
  {
    public:
      DistNode () { dSqr = FLT_MAX; epoch = -1; }
      DistNode (int x, int y, float dSqr, int epoch) 
          : p(x,y) { this->dSqr = dSqr; this->epoch = epoch; }
      bool valid () { return dSqr != FLT_MAX; }
      Point p;
      float dSqr;
      int epoch;
  };

  // compute the distance from every point in the image to the boundary.
  void Segmentation::computeDistanceMap (Image& distanceMap) const
  {
      const int width = getWidth();
      const int height = getHeight();
      distanceMap.resize(width,height);

      Util::Array2D<bool> bmap;
      computeBoundaryMap(bmap);

      Util::Array1D<Point>* frontier = 
          new Util::Array1D<Point> (width*height);
      Util::Array1D<Point>* newFrontier = 
          new Util::Array1D<Point> (width*height);

      int epoch = 0;
      // initialize the frontier set
      Util::Array2D<DistNode> dmap (width,height);
      int frontierSize = 0;
      for (int x = 0; x < width; x++) 
      {
        for (int y = 0; y < height; y++) 
        {
            if (bmap(x,y) == true) 
            {
                (*frontier)(frontierSize++) = Point(x,y);
                dmap(x,y) = DistNode(x,y,0,epoch);
            }
        }
      }

      // enlarge frontier until all pixels are covered
      while (frontierSize > 0) {
          //std::cerr << "frontier size = " << frontierSize << std::endl;
          epoch++;
          int newFrontierSize = 0;
          for (int i = 0; i < frontierSize; i++) {
              const int x = (*frontier)(i).x;
              const int y = (*frontier)(i).y;
              assert (dmap(x,y).valid());
              assert (epoch==1 || bmap(x,y) != true);
              assert (dmap(x,y).epoch >= epoch-1);
              const int hx = dmap(x,y).p.x;
              const int hy = dmap(x,y).p.y;
              assert (dmap(hx,hy).valid());
              assert (bmap(hx,hy) == true);
              for (int v = -1; v <= 1; v++) {
                  for (int u = -1; u <= 1; u++) {
                      const int ix = x + u;
                      const int iy = y + v;
                      if (ix < 0 || ix >= width) { continue; }
                      if (iy < 0 || iy >= height) { continue; }
                      if (u == 0 && v == 0) { continue; }
                      const int dx = hx - ix;
                      const int dy = hy - iy;
                      const float dSqr = dx*dx + dy*dy;
                      if (dSqr < dmap(ix,iy).dSqr) {
                          if (dmap(ix,iy).epoch < epoch) {
                              (*newFrontier)(newFrontierSize++) = 
                                  Point(ix,iy);
                          }
                          dmap(ix,iy) = DistNode(hx,hy,dSqr,epoch);
                      }
                  }
              }
          }
          frontierSize = newFrontierSize;
          Util::swap(frontier,newFrontier);
      }

      // record distance map and homes
      for (int x = 0; x < width; x++) {
          for (int y = 0; y < height; y++) {
              const float d = sqrt(dmap(x,y).dSqr);
              distanceMap(x,y) = d;
          }
      }

      // clean up
      delete frontier;
      delete newFrontier;
  }

  // compute orientation of maximum variance of a list of points, [0,PI)
  float
  primaryOrient (const Array2D<float> points, const int numPoints) 
  {
      assert (numPoints > 0);
      assert (points.size(1) == 2);

      // compute means
      double meanx=0, meany=0;
      for (int i = 0; i < numPoints; i++) 
      {
          const double x = points(i,0);
          const double y = points(i,1);
          meanx += x;
          meany += y;
      }    
      meanx /= numPoints;
      meany /= numPoints;
     
      // compute (unnormalized) covariance matrix
      //	a b
      //  c d
      double a=0, b=0, c=0, d=0;
      for (int i = 0; i < numPoints; i++) 
      {
          const double x = points(i,0);
          const double y = points(i,1);
          a += (x - meanx) * (x - meanx);
          b += (x - meanx) * (y - meany);
          c += (y - meany) * (x - meanx);
          d += (y - meany) * (y - meany);
      }

      // compute eigenvalues
      const double k = sqrt (a*a + 4*b*c - 2*a*d + d*d);
      const double e0 = (a + d - k) / 2;
      const double e1 = (a + d + k) / 2;

      // compute eigenvectors
      const double v0[2] = { -b , a - e0 };
      const double v1[2] = { -b , a - e1 };

      // validate eigenpairs
      const double res1 = a*v0[0] + b*v0[1] - e0*v0[0];
      const double res2 = c*v0[0] + d*v0[1] - e0*v0[1];
      const double res3 = a*v1[0] + b*v1[1] - e1*v1[0];
      const double res4 = c*v1[0] + d*v1[1] - e1*v1[1];
      assert (fabs(res1) < 1e-9);
      assert (fabs(res2) < 1e-9);
      assert (fabs(res3) < 1e-9);
      assert (fabs(res4) < 1e-9);

      // get eigenvector corresponding to largest eigenvalue
      double v[2];
      if (fabs(e0) > fabs(e1)) {
          v[0] = v0[0];
          v[1] = v0[1];
      } else {
          v[0] = v1[0];
          v[1] = v1[1];
      }

      // compute direction of eigenvector
      float theta = atan2(v[1],v[0]);
      assert(finite(theta));
      if (theta < 0) { theta += M_PI; }
      if (theta < 0) { theta = 0; }
      if (theta >= M_PI) { theta = 0; }
      assert (finite(theta));
      assert (theta >= 0 && theta < M_PI);
      return theta;
  }

  // given a boundary map (bmap), estimate the boundary orientation
  // at each boundary point using a window of the given radius
  void Segmentation::computeOrientationMap (const float radius, Image& orientationMap) const
  {
      const int width = getWidth();
      const int height = getHeight();
      const int winrad = (int) ceil (radius);
      Util::Array2D<float> points(width*height,2);

      Util::Array2D<bool> bmap;
      computeBoundaryMap(bmap);

      orientationMap.resize(width,height);
      for (int x = 0; x < width; x++) {
          for (int y = 0; y < height; y++) {
              if (bmap(x,y) != true) { continue; }
              int numPoints = 0;
              for (int v = -winrad; v <= winrad; v++) {
                  for (int u = -winrad; u <= winrad; u++) {
                      int ix = x + u;
                      int iy = y + v;
                      if (ix < 0 || ix >= width) { continue; }
                      if (iy < 0 || iy >= height) { continue; }
                      const float dSqr = u*u + v*v;
                      if (dSqr > radius*radius) { continue; }
                      if (bmap(ix,iy) != true) { continue; }
                      points(numPoints,0) = ix;
                      points(numPoints,1) = iy;
                      numPoints++;
                  }
              }
              const float theta = primaryOrient(points,numPoints);
              orientationMap(x,y) = theta;
          }
      }
  }

  /////////////////////////////////////////////////////////////////////////

  int 
  Segmentation::getWidth () const
  {
      return _map.size(0);
  }

  int 
  Segmentation::getHeight () const
  {
      return _map.size(1);
  }

  int 
  Segmentation::getNumSegments() const
  {
      return _numSegments;
  }

  //
  // get a list of the segment sizes
  //
  void
  Segmentation::segSizes(Array1D<int>& sizes)
  {
    _renumber();
    sizes.resize(_numSegments);
    sizes.init(0);
    for (int x = 0; x < _map.size(0); x++) 
    {
      for (int y = 0; y < _map.size(1); y++) 
      {
        sizes(_map(x,y))++;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////

  const char*
  Segmentation::getDate () const
  {
      return _date;
  }

  const char* 
  Segmentation::getUser () const
  {
      return _user;
  }

  const char* 
  Segmentation::getImage () const
  {
      return _image;
  }

  int 
  Segmentation::getImageID () const
  {
      return _imageID;
  }

  int 
  Segmentation::getUserID () const
  {
      return _userID;
  }

  int 
  Segmentation::getPresentation () const
  {
      return _presentation;
  }
  bool 
  Segmentation::isColor () const
  {
      return ! isGray ();
  }

  bool 
  Segmentation::isGray () const
  {
      return _presentation & presGray;
  }

  bool 
  Segmentation::isInverted () const
  {
      return _presentation & presInvert;
  }

  bool 
  Segmentation::isFlipFloped () const
  {
      return _presentation & presFlipFlop;
  }

}
