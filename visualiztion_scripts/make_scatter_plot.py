# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 13:52:36 2012
Plots distribution of desired isotope (Iso) for all times in a pfltoran h5 file.(inputfile);
@author: wpgardn
"""

from pflotran_tools import *
from pylab import *
import h5py
import pdb

# the name of the input file is needed below - pflotran.h5 is the default
inputfile= str(raw_input('Please enter the name of output file: '));
i = 0;
s_bool = raw_input('Save all figures (Y/N): ');
c_flag = 'Y';
dd = ReadPflotranObsPt(inputfile);
for k in dd.dtype.names:
    print k;
Iso = str(raw_input('Type in variable to plot from list above (must match exactly): '));

#pdb.set_trace();
figure();
while c_flag == 'Y':
    loglog(dd[dd.dtype.names[0]],dd[Iso],label=Iso);
    xlabel('X (m)')
    ylabel('Variable')
    legend(loc='best');
    c_flag = raw_input('Thank you Sir/Madame, Would you like to plot another? (Y/N): ');
    if c_flag == 'Y':
       for k in dd.dtype.names:
          print k;
       Iso = str(raw_input('Type in variable to plot from list above (must match exactly): ')); 
    if s_bool == 'Y':
        savefig('scatter.png');

show();
