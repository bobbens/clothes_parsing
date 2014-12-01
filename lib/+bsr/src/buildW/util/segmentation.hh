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

#ifndef SEGMENTATION_HH 
#define SEGMENTAITON_HH 

#include <stdio.h>
#include <assert.h>

#include "image.hh"

namespace Util
{
  class Segmentation 
  {
    public:
      enum FileFormat 
      {
        ffASCII	= 0x1,	// ASCII / binary
        ffCr	= 0x2	// cr / map
      };

      enum Presentation 
      {
        presGray 	= 0x1,
        presInvert 	= 0x2,
        presFlipFlop 	= 0x4
      };

      static const int defaultFormat = ffASCII | ffCr;

      Segmentation();
      Segmentation(const Array2D<int> parititon);
      Segmentation(const int width, const int height);
      ~Segmentation();

      void readFile(const char* file);
      void writeFile(const char* file, int format = defaultFormat) const;

      //get a list of segment sizes
      void segSizes(Array1D<int>& sizes);

      //merge segments smaller than dustthresh into nearby segments
      void mergeDust(int dustthresh);

      //split disconnected segments into unique segments.
      void fragment(float radius = 1.5);

      //create a boundary map from the segmentation 
      void computeBoundaryMap(Util::Array2D<bool>& boundaryMap) const;
      void computeBoundaryMapHalf(Util::Array2D<bool>& boundaryMap) const;


      //compute distance from each pixel to nearest boundary
      void computeDistanceMap(Image& distanceMap) const;

      //assign an orientation to every boundary point
      void computeOrientationMap(const float radius, Image& orientationMap) const;

      //segmentation dimensions
      int getWidth () const;
      int getHeight () const;
      int getNumSegments() const;

      //Array2D style access
      int& operator()(int x, int y)
      {
        return _map(x,y);
      }
      const int& operator()(int x, int y) const
      {
        return _map(x,y);
      }

      //access segmentation database details
      const char* getDate () const;
      const char* getUser () const;
      const char* getImage () const;
      int getImageID () const;
      int getUserID () const;
      int getPresentation () const;
      bool isColor () const;
      bool isGray () const;
      bool isInverted () const;
      bool isFlipFloped () const;

    protected:

      //renumber segments
      void _renumber();

      //flood fill type operation used in finding connected components
      void _replace (int row, int col, int oldId, int newId, float radius);

      Array2D<int> _map;
      int _numSegments;

      char*	_date;
      char*	_image;
      char*	_user;
      int _imageID;
      int _userID;
      int _presentation;
   


    friend std::istream & operator >> (std::istream & in, Segmentation& s)
    {
      in >> s._map;
      s._date = NULL;
      s._image = NULL;
      s._user = NULL;
      s._imageID = 0;
      s._userID = 0;
      s._presentation = 0;
      s._renumber();
      return in;
    }


    friend std::ostream & operator << (std::ostream & out, const Segmentation& s)
    {
      out << s._map;
      return out;
    }
  }; //class Segmentation
} //namespace Util
#endif // SEGMENTAITON_HH 




