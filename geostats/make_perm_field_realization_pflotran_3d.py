# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 13:53:52 2012
This script will create a hypothetical permeability distribution.  First will create a spatial seed drawn from a 
hypothetical distribution 
@author: wpgardn
"""

#import sys
#sys.path.append(r'../shared')

from numpy import *
from geo_bsd import *
from matplotlib import *
from pylab import *
from scipy import *
from scipy.stats import lognorm
from matplotlib.colors import LogNorm
#from gslib import *
#from grid_3d import *
#from gaussian_cdf import *
#from variogram_routines import *
from geo_bsd.routines import *
import pdb
import h5py

# number of cells
i_max = 500
j_max = 500
k_max = 50

#grid dimensions
xdim = 500.;
ydim = 500.;
zdim = 50.;

# create the empty data property
prop_ijk = zeros((i_max, j_max, k_max));
prop_ijk = require(prop_ijk, dtype=float32, requirements=['F']);

array_defined = zeros((i_max, j_max, k_max));
array_defined = require(array_defined, dtype=uint8, requirements=['F']);

initial_data = copy(prop_ijk);

#Now create an emtpy property
prop_initial = (initial_data, array_defined);
prop_initial = ContProperty(prop_initial[0], prop_initial[1]);

# Define 3D grid 
grid = SugarboxGrid(i_max, j_max, k_max)

#lets try making a hypothetical sampling of data with no-spatial information
#one order of magnitude
mu2=-40.
sigma2 = 10.1;
#two orders of magnitude?
#mu2=-50.
#sigma2=20.1

numargs = lognorm.numargs
[s] = [0.1]*numargs
#this will grab 400 samples from the distribution with paramters given above - see perm_dist.py to
#see what this distribution looks like.
rv = lognorm(s,loc=mu2,scale=sigma2);
sts = rv.stats()
print 'mean is = ' + str(exp(sts[0]));
print 'std is = ' + str(sqrt(sts[1]));
prop_init_forcdf = lognorm.rvs(s,loc=mu2,scale=sigma2,size=2000);
a = exp(reshape(prop_init_forcdf,(20,20,5)));
b = ones((20,20,5));
figure();
subplot(1,3,2);
title('Seed Data');
bins = logspace(-15,-8,num=150);
hist(exp(prop_init_forcdf),bins=bins);
gca().set_xscale("log");

#since I'm not try to condition the data I think I can just add in the data at any point on the grid
prop_initial.data[:20,:20,:5]+=a;
prop_initial.mask[:20,:20,:5]+=b;



#create a varigram object
#correlation lengths in meters;
ycorrlen=20;
xcorrlen=20;
zcorrlen=20;
#correlation lengths in grid dimensions:
ycorr=int(ycorrlen*(i_max/ydim));
xcorr=int(xcorrlen*(j_max/xdim));
zcorr=int(zcorrlen*(k_max/zdim));
#search radii in grid dimensions:
yrad = ycorr*2;
xrad = xcorr*2;
zrad = zcorr*2;

variogram2 = CovarianceModel(type=covariance.spherical, ranges = (ycorr, xcorr, zcorr), sill = 1, angles = (0, 0, 0))
sgs_params = {"cov_model": variogram2, "radiuses": (yrad, xrad, zrad), "max_neighbours": 20 }

#test to see if we can do unconditional sgs now...
sgs_initial_data_uncond = sgs_simulation(prop_initial, grid, use_harddata=False, cdf_data=calc_cdf(prop_initial), seed=54783, **sgs_params)

# Simulated conditional realization in initial
figure()
#subplot(1,3,2);
#pdb.set_trace();
imshow(sgs_initial_data_uncond[0][:,:,2], norm=LogNorm(vmin = a.min(), vmax = a.max()));
t = logspace(log10(a.min()),log10(a.max()),num=10);
bar = colorbar(ticks=t,format='%.2e');
bar.set_label('Permeability Field');
title("Simulated unconditional realization")


figure()
#subplot(1,3,3);
title('Simulated Data');
bins = logspace(-15,-8,num=150);
hist(ravel(sgs_initial_data_uncond[0][:,:,:]),bins=bins);
gca().set_xscale("log");


#now need to figure out how to write sgs_initial_data to hdf5 in right format.
# pflotran.h5 is a 1d array 0 is lower left and numbering increases by right hand rule: i.e. (x,y,z) i,j,k
# index = i + nx*j + nx*ny*k
# so
#sgs_realization_pflotran = ndarray([]);
s_flag=1
for k in xrange(k_max): 
#    pdb.set_trace();
    for j in xrange(j_max-1,-1,-1): # (0,0) (x,y) is top left so we want to start with last row and work up
        if s_flag == 1:
            sgs_realization_pflotran = sgs_initial_data_uncond[0][j,:,k];
            s_flag=0;
        else:
            sgs_realization_pflotran = concatenate((sgs_realization_pflotran,sgs_initial_data_uncond[0][j,:,k]));

# now create the h5 array - this is a little klugey, but it works i think
f = h5py.File('perm_field.h5','w');
#f = h5py.File('por_field.h5','w');
dset = f.create_dataset('Permeability',data=sgs_realization_pflotran);
#dset = f.create_dataset('Porosity',data=sgs_realization_pflotran);
f['Cell Ids']=range(1,len(sgs_realization_pflotran)+1);
#f['Porosity'][:]=f['Porosity'][:]/100.;
f.close();
show();
