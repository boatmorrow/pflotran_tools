# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 10:40:08 2013

@author: wpgardn
"""

from pflotran_2d_vizualization import *
from pylab import *
import h5py
import pdb
import cfc_tools as cfc
import sf6_tools as sf6


#input file
ifile = 'pflotran.h5'

#model year start
yr_start = 1930;
output_yr = 2000;
x_dist = 300;

#get an age slice in 2000 at x=250 m 
sxz,yy,zz = GetTracerYZSlice('Total Tracer_Age [M]',x_dist,output_yr-yr_start,ifile);
sxz = sxz/3.156e7
print 'Average Modeled Mean Age = '+str(mean(sxz)) + ' years'
print 'Standard Deviation of Mean Age is '+str(std(sxz))+ ' years'

#get an 3H slice
H3_sxz,yy,zz = GetTracerYZSlice('Total H3 [M]',x_dist,output_yr-yr_start,ifile);

#get an 3He slice 
He3_sxz,yy,zz = GetTracerYZSlice('Total He3 [M]',x_dist,output_yr-yr_start,ifile);

#get an CFC slice
CFC12_sxz,yy,zz = GetTracerYZSlice('Total CFC12 [M]',x_dist,output_yr-yr_start,ifile);

#get an SF6 slice
SF6_sxz,yy,zz = GetTracerYZSlice('Total SF6 [M]',x_dist,output_yr-yr_start,ifile);

#get and Ar39 slice
Ar39_sxz,yy,zz = GetTracerYZSlice('Total Ar39 [M]',x_dist,output_yr-yr_start,ifile);

#calculate 3H/3He piston flow age
lambda_h3=0.0558;
c_he_atm = 2.8253e-15;
pist_age = 1./lambda_h3*log((He3_sxz-c_he_atm)/H3_sxz+1.);
pist_age[where(pist_age>70)]=70;
print 'Average 3H/3He Piston Flow Age = '+str(mean(pist_age)) + ' years'
print 'Standard Deviation of 3H/3He Piston Flow Age is '+str(std(pist_age))+ ' years'

#calculate the CFC12 piston flow age

#calculate CFC12 age
vfunc = vectorize(cfc.equil_air_conc);
vfunc2 = vectorize(cfc.age_date);
eac = vfunc(CFC12_sxz,25.0,0.0)/1e-12; #pptv
pist_year_cfc12 = vfunc2(eac);
pist_age_cfc12 = output_yr-pist_year_cfc12;
pist_age_cfc12[where(pist_age_cfc12<0)]=0;
print 'Average CFC12 Piston Flow Age = '+str(mean(pist_age_cfc12)) + ' years'
print 'Standard Deviation of CFC12 Piston Flow Age is '+str(std(pist_age_cfc12))+ ' years'

#calculate the SF6 age
vfunc = vectorize(sf6.equil_air_conc);
vfunc2 = vectorize(sf6.age_date);
eac = vfunc(SF6_sxz,25.0,0.0)/1e-15; #pptv
pist_year_sf6 = vfunc2(eac);
pist_age_sf6 = output_yr-pist_year_sf6;
pist_age_sf6[where(pist_age_sf6<0)]=0;
print 'Average SF6 Piston Flow Age = '+str(mean(pist_age_sf6)) + ' years'
print 'Standard Deviation of SF6 Piston Flow Age is '+str(std(pist_age_sf6))+ ' years'

#calculate Ar39 piston flow age
lambda_Ar39=2.58e-3;
c_Ar39_atm = 1.0211e-20;
pist_age_Ar39 = -lambda_Ar39**-1*log(Ar39_sxz/c_Ar39_atm);
#pdb.set_trace();
pist_age_Ar39[where(pist_age_Ar39<100)]=0.01;
#pist_age_Ar39[where(pist_age_Ar39>1000)]=1000;
print 'Average Ar39 Piston Flow Age = '+str(mean(pist_age_Ar39)) + ' years'
print 'Standard Deviation of Ar39 Piston Flow Age is '+str(std(pist_age_Ar39))+ ' years'

#pdb.set_trace();
figure();
#plot(yy,(sxz/3.156e7)[:,0]);
subplot(5,1,1);
hist((sxz),bins=50,label='Mean Age',color='g');
xmin,xmax = xlim();
xmin=0
xlim(xmin,xmax);
legend();
subplot(5,1,2);
hist(pist_age,bins=1,label=r'$^3$H/$^3$He Piston Flow Age',color='r');
xlim(xmin,xmax);
legend();
subplot(5,1,3);
hist(pist_age_cfc12,bins=1,label=r'CFC12 Piston Flow Age',color='b');
xlim(xmin,xmax);
legend();
subplot(5,1,4);
hist(pist_age_sf6,bins=1,label=r'SF6 Piston Flow Age',color='k');
xlim(xmin,xmax);
legend();
subplot(5,1,5);
hist(pist_age_Ar39,bins=50,label=r'Ar39 Piston Flow Age',color='m');
#pdb.set_trace();
xlim(xmin,xmax);
xlabel(r'$\bar{A} (yr)$')
legend();
show();