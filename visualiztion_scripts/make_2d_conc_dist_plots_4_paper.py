# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 15:40:57 2013

@author: wpgardn
"""
from matplotlib import *
from pflotran_2d_vizualization import *
import pdb

ifile = 'pflotran.h5';


rcParams['figure.figsize']=(14,10);

yr_start = 1930;
output_yr = 2000;

figure();
subplot(2,2,1);
#PlotMeanAgeDistribution(output_yr-yr_start,ifile);
PlotTracerDistribution('Total H3 [M]',output_yr-yr_start,ifile);

subplot(2,2,2);
#Plot3He3HAgeDistriubtion(output_yr-yr_start,ifile);
PlotTracerDistribution('Total He3 [M]',output_yr-yr_start,ifile);

subplot(2,2,3);
#PlotCFC12AgeDistriubtion(output_yr-yr_start,ifile,t_start=yr_start);
PlotTracerDistribution('Total CFC12 [M]',output_yr-yr_start,ifile);

subplot(2,2,4);
#PlotSF6AgeDistriubtion(output_yr-yr_start,ifile,t_start=yr_start);
PlotTracerDistribution('Total SF6 [M]',output_yr-yr_start,ifile);

figure();
subplot(2,2,1);
#PlotMeanAgeDistribution(output_yr-yr_start,ifile);
PlotTracerDistribution('Total Ar39 [M]',output_yr-yr_start,ifile);

subplot(2,2,2);
#Plot3He3HAgeDistriubtion(output_yr-yr_start,ifile);
PlotTracerDistribution('Total Kr81 [M]',output_yr-yr_start,ifile);

subplot(2,2,3);
#PlotCFC12AgeDistriubtion(output_yr-yr_start,ifile,t_start=yr_start);
PlotTracerDistribution('Total He4 [M]',output_yr-yr_start,ifile);

#subplot(2,2,4);
#PlotSF6AgeDistriubtion(output_yr-yr_start,ifile,t_start=yr_start);
#PlotTracerDistribution('Total SF6 [M]',output_yr-yr_start,ifile);

show();
