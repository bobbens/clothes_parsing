#!/usr/bin/env python2
from nltk.corpus import wordnet as wn
import numpy as np
import scipy.io

synset_names = [
   'accessory.n.01',
   'bag.n.04',
   'belt.n.02',
   'blazer.n.01',
   'blouse.n.01',
   'leotard.n.01', # BODYSUIT, NOT IN FASHIONISTA 0.3
   'boot.n.01',
   'bra.n.01',
   'bracelet.n.01',
   'cape.n.02',
   'cardigan.n.01',
   'clog.n.01',
   'coat.n.01',
   'dress.n.01',
   'earring.n.01',
   'flats.n.01',
   'glasses.n.01',
   'glove.n.02',
   'hair.n.01',
   'hat.n.01',
   'platform.n.05', # should be heel
   'jacket.n.01',
   'underwear.n.01', # INTIMATE, NOT IN FASHIONISTA 0.3
   'jean.n.01',
   'jumper.n.03',
   'legging.n.01',
   'loafer.n.02',
   'necklace.n.01',
   'pants.n.01',
   'pump.n.03',
   'purse.n.01',
   'ring.n.08',
   'romper.n.02',
   'sandal.n.01',
   'scarf.n.01',
   'shirt.n.01',
   'shoe.n.01',
   'shorts.n.01',
   'skin.n.01',
   'skirt.n.01',
   'sneaker.n.01',
   'sock.n.01',
   'stocking.n.01',
   'suit.n.01',
   'sunglasses.n.01',
   'sweater.n.01',
   'sweatshirt.n.01',
   't-shirt.n.01',
   'tie.n.01',
   'tights.n.01',
   'top.n.10',
   'vest.n.01',
   'wallet.n.01',
   'wristwatch.n.01',
   'wedgie.n.01'
   ]

synset_list = []
for i in range(len(synset_names)):
   ss = wn.synset( synset_names[i] )
   synset_list.append( ss )
   #hyp = lambda s:s.hypernyms()
   #from pprint import pprint
   #pprint(ss.tree(hyp,5))



#ic    = wn.ic()
n        = len(synset_list)
pathmat  = np.zeros((n,n))
wupmat   = np.zeros((n,n))
lchmat   = np.zeros((n,n))
dmat     = np.zeros((n,n))
for i1, s1 in enumerate(synset_list):
   for i2, s2 in enumerate(synset_list):
      pathmat[i1,i2] = s1.path_similarity(s2)
      wupmat[i1,i2]  = s1.wup_similarity(s2)
      lchmat[i1,i2]  = s1.lch_similarity(s2)
      dmat[i1,i2]    = s1.shortest_path_distance(s2)
      #dmat[i1,i2] = jcn_similarity( s1, s2, ic )

scipy.io.savemat( 'wn_pathmat_v2.mat', mdict={'wn_sim':pathmat} )
scipy.io.savemat( 'wn_wupmat_v2.mat', mdict={'wn_sim':wupmat} )
scipy.io.savemat( 'wn_lchmat_v2.mat', mdict={'wn_sim':lchmat} )
scipy.io.savemat( 'wn_distmat_v2.mat', mdict={'wn_dist':dmat} )
np.set_printoptions(threshold='nan')
#print(dmat)
#print(n)

