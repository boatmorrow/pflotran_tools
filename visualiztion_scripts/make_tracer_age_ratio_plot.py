# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 10:03:17 2013
plot tracer age/mean age ratio vs. unitless distance
@author: wpgardn
"""

from pflotran_2d_vizualization import *
from pylab import *
import h5py
import pdb
import cfc_tools as cfc
import sf6_tools as sf6
import matplotlib as mpl

mpl.rcParams['axes.labelsize'] = 24;
mpl.rcParams['figure.figsize']=(14,10)

#input file
ifile = 'pflotran.h5'

#model year start
yr_start = 1930;
output_yr = 2000;
x_dist_vec = [3.,5.,7.,10.,15.,20.,50.,75.,100.,200.,300.,400.,500.];
MM_age = [];
MM_std = [];
H3He3M_age = [];
H3He3_std = [];
CFC12M_age = [];
CFC12_std = [];
SF6_age = [];
SF6_std = [];
Ar39_age = [];
Ar39_std = [];

#get an age slice in 2000 at x=250 m 
for x_dist in x_dist_vec:
    sxz,yy,zz = GetTracerYZSlice('Total Tracer_Age [M]',x_dist,output_yr-yr_start,ifile);   
    sxz = sxz/3.156e7
    MM_age.append(mean(sxz));
    MM_std.append(std(sxz));
    #print 'Average Modeled Mean Age = '+str(mean(sxz)) + ' years'
    #print 'Standard Deviation of Mean Age is '+str(std(sxz))+ ' years'

    #get an 3H slice
    H3_sxz,yy,zz = GetTracerYZSlice('Total H3 [M]',x_dist,output_yr-yr_start,ifile);
    
    #get an 3He slice 
    He3_sxz,yy,zz = GetTracerYZSlice('Total He3 [M]',x_dist,output_yr-yr_start,ifile);

    #get an CFC slice
    CFC12_sxz,yy,zz = GetTracerYZSlice('Total CFC12 [M]',x_dist,output_yr-yr_start,ifile);

    #get an SF6 slice
    SF6_sxz,yy,zz = GetTracerYZSlice('Total SF6 [M]',x_dist,output_yr-yr_start,ifile);
    
    #get an Ar39 slice
    Ar39_sxz,yy,zz = GetTracerYZSlice('Total Ar39 [M]',x_dist,output_yr-yr_start,ifile);

    #calculate 3H/3He piston flow age
    lambda_h3=0.0558;
    c_he_atm = 2.8253e-15;
    pist_age = 1./lambda_h3*log((He3_sxz-c_he_atm)/H3_sxz+1.);
    pist_age[where(pist_age>70)]=70;
    H3He3M_age.append(mean(pist_age));
    H3He3_std.append(std(pist_age));
    #print 'Average 3H/3He Piston Flow Age = '+str(mean(pist_age)) + ' years'
    #print 'Standard Deviation of 3H/3He Piston Flow Age is '+str(std(pist_age))+ ' years'

    #calculate the CFC12 piston flow age

    #calculate CFC12 age
    vfunc = vectorize(cfc.equil_air_conc);
    vfunc2 = vectorize(cfc.age_date);
    eac = vfunc(CFC12_sxz,25.0,0.0)/1e-12; #pptv
    pist_year_cfc12 = vfunc2(eac);
    pist_age_cfc12 = output_yr-pist_year_cfc12;
    pist_age_cfc12[where(pist_age_cfc12<0)]=0;
    CFC12M_age.append(mean(pist_age_cfc12));
    CFC12_std.append(std(pist_age_cfc12));
    #print 'Average CFC12 Piston Flow Age = '+str(mean(pist_age_cfc12)) + ' years'
    #print 'Standard Deviation of CFC12 Piston Flow Age is '+str(std(pist_age_cfc12))+ ' years'
    
    #calculate the SF6 age
    vfunc = vectorize(sf6.equil_air_conc);
    vfunc2 = vectorize(sf6.age_date);
    eac = vfunc(SF6_sxz,25.0,0.0)/1e-15; #pptv
    pist_year_sf6 = vfunc2(eac);
    pist_age_sf6 = output_yr-pist_year_sf6;
    pist_age_sf6[where(pist_age_sf6<0)]=0;
    SF6_age.append(mean(pist_age_sf6));
    SF6_std.append(std(pist_age_sf6));
    #print 'Average SF6 Piston Flow Age = '+str(mean(pist_age_sf6)) + ' years'
    #print 'Standard Deviation of SF6 Piston Flow Age is '+str(std(pist_age_sf6))+ ' years'
    
    #calculate Ar39 piston flow age
    lambda_h3=0.0558;
    c_he_atm = 2.8253e-15;
    lambda_Ar39=2.58e-3;
    c_Ar39_atm = 1.0211e-20;
    pist_age_Ar39 = -lambda_Ar39**-1*log(Ar39_sxz/c_Ar39_atm);
    #pdb.set_trace();
    pist_age_Ar39[where(pist_age_Ar39<100)]=0;
    pist_age_Ar39[where(pist_age_Ar39>1000)]=1000;
    Ar39_age.append(mean(pist_age_Ar39));
    Ar39_std.append(std(pist_age_Ar39));
    #print 'Average 3H/3He Piston Flow Age = '+str(mean(pist_age)) + ' years'
    #print 'Standard Deviation of 3H/3He Piston Flow Age is '+str(std(pist_age))+ ' years'
    

MM_age = array(MM_age);
MM_std = array(MM_std);
H3He3M_age = array(H3He3M_age);
H3He3_std = array(H3He3_std);
CFC12M_age = array(CFC12M_age);
CFC12_std = array(CFC12_std);
SF6_age = array(SF6_age);
SF6_std = array(SF6_std);
Ar39_age = array(Ar39_age);
Ar39_std = array(Ar39_std);

CFC12_max_age = 60.;
H3He3_max_age = 70.;
SF6_max_age = 45.;
Ar39_max_age = 1000.;
vx_sxz,yy,zz = GetTracerYZSlice('Liquid X-Velocity',250.,output_yr-yr_start,ifile);
mean_vx = mean(vx_sxz);
print 'Average Velocity = ' +str(mean_vx)+ ' m/yr';

fig1=figure();
#pdb.set_trace();
subplot(2,1,1);
plot(x_dist_vec/(mean_vx*H3He3_max_age),H3He3M_age/MM_age,'ro',label=r'$^3$H/$^3$He',ms=8);
plot(x_dist_vec/(mean_vx*CFC12_max_age),CFC12M_age/MM_age,'bs',label=r'CFC12',ms=8);
plot(x_dist_vec/(mean_vx*SF6_max_age),SF6_age/MM_age,'k^',label=r'SF$_6$',ms=8);
plot(x_dist_vec/(mean_vx*Ar39_max_age),Ar39_age/MM_age,'mD',label=r'$^{39}Ar$',ms=8);
#xlabel(r'$\huge{\frac{X}{v_{x_{av}}\dot t_{tr_{max}}}}$');
xlabel(r'$\frac{X}{v_{x_{av}} \cdot t_{tr_{max}}}$');
ylabel(r'$\bar{A_{tr}} / \bar{A_{m}}$')
legend(loc='best',numpoints=1);

subplot(2,1,2);
plot(x_dist_vec,H3He3M_age/MM_age,'ro',label=r'$^3$H/$^3$He',ms=8);
plot(x_dist_vec,CFC12M_age/MM_age,'bs',label=r'CFC12',ms=8);
plot(x_dist_vec,SF6_age/MM_age,'k^',label=r'SF$_6$',ms=8);
plot(x_dist_vec,Ar39_age/MM_age,'mD',label=r'$^{39}Ar$',ms=8);
xlabel(r'$X$');
ylabel(r'$\bar{A_{tr}} / \bar{A_{m}}$')
#legend(loc='best',numpoints=1);
fig1.subplots_adjust(hspace=.3);

fig2=figure();
#pdb.set_trace();
subplot(2,1,1);
plot(x_dist_vec/(mean_vx*H3He3_max_age),H3He3_std/MM_std,'ro',label=r'$^3$H/$^3$He',ms=8);
plot(x_dist_vec/(mean_vx*CFC12_max_age),CFC12_std/MM_std,'bs',label=r'CFC12',ms=8);
plot(x_dist_vec/(mean_vx*SF6_max_age),SF6_std/MM_std,'k^',label=r'SF$_6$',ms=8);
plot(x_dist_vec/(mean_vx*Ar39_max_age),Ar39_std/MM_std,'mD',label=r'$^{39}Ar$',ms=8);
#xlabel(r'$\huge{\frac{X}{v_{x_{av}}\dot t_{tr_{max}}}}$');
xlabel(r'$\frac{X}{v_{x_{av}} \cdot t_{tr_{max}}}$');
ylabel(r'$\bar{\sigma_{tr}} / \bar{\sigma_{m}}$');
ylim(-.03,1.5);
legend(loc='best',numpoints=1);

subplot(2,1,2);
plot(x_dist_vec,H3He3_std/MM_std,'ro',label=r'$^3$H/$^3$He',ms=8);
plot(x_dist_vec,CFC12_std/MM_std,'bs',label=r'CFC12',ms=8);
plot(x_dist_vec,SF6_std/MM_std,'k^',label=r'SF$_6$',ms=8);
plot(x_dist_vec,Ar39_std/MM_std,'mD',label=r'$^{39}Ar$',ms=8);
xlabel(r'$X$');
ylabel(r'$\bar{\sigma_{tr}} / \bar{\sigma_{m}}$');
ylim(-.05,1.5);
#legend(loc='best',numpoints=1);
fig2.subplots_adjust(hspace=.3);


#xlim(0,5);
show();

