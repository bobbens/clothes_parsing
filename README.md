Clothes Parsing
===============

Overview
--------

This code provides an implementation of the research paper:

```
  A High Performance CRF Model for Clothes Parsing
  Edgar Simo-Serra, Sanja Fidler, Francesc Moreno-Noguer, and Raquel Urtasun
  Asian Conference on Computer Vision (ACCV), 2014
```

The code here allows training and testing of a model that got state-of-the-art
results on the
[Fashionista](http://vision.is.tohoku.ac.jp/~kyamagu/ja/research/clothing_parsing/)
dataset at the time of publication.


License
-------

```
  Copyright (C) <2014> <Edgar Simo-Serra, Sanja Fidler, Francesc Moreno-Noguer, Raquel Urtasun>

  This work is licensed under the Creative Commons
  Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy
  of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or
  send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

  Edgar Simo-Serra, Institut de Robotica i Informatica Industrial (CSIC/UPC), December 2014.
  esimo@iri.upc.edu, http://www-iri.upc.es/people/esimo/
```


Installation
------------

In order to get started first checkout out the source code and then extract the
features:

```
# Check out the git and cd into it as working directory
git clone https://github.com/bobbens/clothes_parsing.git
cd clothes_parsing
# Get and unpack the necessary features
wget http://www.iri.upc.edu/people/esimo<%= item_named('/data/poseseg/').path%>
tar xvjf poseseg.tar.bz2 
```

The
[dSP](http://www.alexander-schwing.de/projectsGeneralStructuredPredictionLatentVariables.php)
dependency must also be compiled. This can be done by:

```
cd lib/dSP_5.1
make # First edit the Makefile if necessary
```

Usage
-----

You can reproduce results simply by running from Matlab:

```
sm = segmodel( 'PROFILE', '0.16', 'use_real_pose', false ); % Load the model, parameters can be set here
sm = sm.train_misc_unaries(); % Trains some misc stuff
sm = sm.train_MRF(); % Actually sets up and trains the CRF
R = sm.test_MRF_segmentation() % Performs testing and outputs results
```

This should generate an output like:

```
 BUILDING MRF OUTPUT 29 CLASSES (REAL POSE=0)...
 UNARIES:
    bgbias
    logreg:       29
    cpmc_logreg:  29
    cpmc
    shapelets
 HIGHER ORDER
    similarity
    limbs
 Initializing Image 011 / 350...   0.4 seconds!   

 ...

 Tested MRF in 319.0 seconds
 350 / 350... 

 R = 

     confusion: [29x29 double]
     order: [29x1 double]
     acc: 0.8432
     pre: [29x1 double]
     rec: [29x1 double]
     f1: [29x1 double]
     voc: [29x1 double]
     avr_pre: 0.3007
     avr_rec: 0.3292
     avr_f1: 0.3039
     avr_voc: 0.2013
```

Please note that due to stochastic components and differences between software
versions, the numbers will not be exactly the same as the paper. For the paper
all results were obtained on a linux machine running Ubuntu 12.04 with Matlab
R2012a (7.14.0.739) 64-bit (glnxa64).

If you use this code please cite:

```
 @InProceedings{SimoSerraACCV2014,
    author = {Edgar Simo-Serra and Sanja Fidler and Francesc Moreno-Noguer and Raquel Urtasun},
    title = {{A High Performance CRF Model for Clothes Parsing}},
    booktitle = "Proceedings of the Asian Conference on Computer Vision (2014)",
    year = 2014
 }
```

Acknowledgments
---------------

We would like to give our thanks to [Kota
Yamaguchi](http://vision.is.tohoku.ac.jp/~kyamagu/) for his excellent code
which we have used as a base for our model.

The different codes we have used (in alphabetical order):

 * [clothing_parsing](http://vision.is.tohoku.ac.jp/~kyamagu/ja/research/clothing_parsing/)
 * [CPMC](http://www.maths.lth.se/matematiklth/personal/sminchis/code/cpmc/index.html)
 * [dSP](http://www.alexander-schwing.de/projectsGeneralStructuredPredictionLatentVariables.php)
 * [gPb](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)
 * [grTheory](http://www.mathworks.com/matlabcentral/fileexchange/4266-grtheory-graph-theory-toolbox)
 * [liblinear](http://www.csie.ntu.edu.tw/~cjlin/liblinear/)
 * [Pose](http://www.ics.uci.edu/~dramanan/software/pose/)
 * [Selective Search](http://koen.me/research/selectivesearch/)


Changelog
---------

December 2014: Initial version 1.0 release




