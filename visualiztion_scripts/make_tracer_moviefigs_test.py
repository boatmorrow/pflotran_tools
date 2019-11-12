# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 13:52:36 2012
Plots distribution of desired isotope (Iso) for all times in a pfltoran h5 file.(inputfile);
@author: wpgardn
"""

from pflotran_2d_vizualization import *
from pylab import *
import h5py
import pdb

# the name of the input file is needed below - pflotran.h5 is the default
inputfile= str(raw_input('Please enter the name of output file '));

#Iso='Liquid Pressure [Pa]'
#Iso='Total I_129(aq) [M]'
#Iso='Total U_233(aq) [M]'
#Iso='Total Am_241(aq) [M]'
# loop through file
# get the grid
xx,yy,zz = GetCellCenters(inputfile);
times,timekeys = GetTimeInfo(inputfile);
f1=h5py.File(inputfile);
i = 0;
#pdb.set_trace();
for k in f1[timekeys[0]].keys():
    print k;
Iso = str(raw_input('Select field to image from list above '));

for k in timekeys:
    figure(figsize=(16,8));
    #imshow(transpose(log10(f1[k]['Total I_129(aq) [M]'][:,0,:])),origin='lower',extent=(0,max(xx),0,max(zz)),aspect='auto'); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    imshow(transpose(f1[k][Iso][:,0,:]),origin='lower',extent=(0,max(xx),0,max(zz)),aspect='auto'); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    #bar = colorbar();
    bar.set_label(Iso);
    xlabel('X (m)')
    ylabel('Z (m)')
    labeltxt='Time =' + str(times[i]) + 'yr'
    title(labeltxt);
    i = i+1;
    show();
