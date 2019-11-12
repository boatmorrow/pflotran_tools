# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 13:52:36 2012
create tracer concentration plots at a point
@author: wpgardn
"""

from pflotran_2d_vizualization import *
from pylab import *
import h5py
import pdb
import cfc_tools as cfc
import sf6_tools as sf6

mpl.rcParams['axes.labelsize'] = 18;
mpl.rcParams['figure.figsize']=(14,10)

# locaction of point
x = 50.;
y = 10.;
z = 0.5;

#inputfile
ifile = 'pflotran.h5'

f1 = h5py.File(ifile,'r');

#for k in f1[f1.keys()[1]].keys() :
#    try:
#        k.split()[2];
#    except IndexError:
#        continue;
#    else:
#        if k.split()[2]=='[M]':
#            t_i,c_i = GetTracerObsSeries(k,x,y,z,ifile); 
#            #pdb.set_trace()
#            semilogy(t_i,c_i,'o',label=k);

#legend();
#show();

# make a plot of modeled age and 3h/3He age
t_obs = 2012
# Get modeled mean age
t_i,A_i = GetTracerObsSeries('Total Tracer_Age [M]',x,y,z,ifile);
A_i = A_i/3.156e7; #years

#get tracer concentrations of interest
t_i,H3 = GetTracerObsSeries('Total H3 [M]',x,y,z,ifile);
t_i,He3 = GetTracerObsSeries('Total He3 [M]',x,y,z,ifile);
t_i,CFC12 = GetTracerObsSeries('Total CFC12 [M]',x,y,z,ifile);
t_i,SF6 = GetTracerObsSeries('Total SF6 [M]',x,y,z,ifile);

#calculate 3H/3He age
t_i = t_i+1930; #years
lambda_h3=0.0558;
c_he_atm = 2.8253e-15;
pist_age = 1./lambda_h3*log((He3-c_he_atm)/H3+1.);

#calculate CFC12 age
vfunc = vectorize(cfc.equil_air_conc);
vfunc2 = vectorize(cfc.age_date);
eac = vfunc(CFC12,25.0,0.0)/1e-12; #pptv
pist_year_cfc12 = vfunc2(eac);
pist_age_cfc12 = t_i-pist_year_cfc12;
pist_age_cfc12[where(pist_age_cfc12<0)]=0;
#calculate the SF6 age
vfunc = vectorize(sf6.equil_air_conc);
vfunc2 = vectorize(sf6.age_date);
eac = vfunc(SF6,25.0,0.0)/1e-15; #pptv
pist_year_sf6 = vfunc2(eac);
pist_age_sf6 = t_i-pist_year_sf6;
pist_age_sf6[where(pist_age_sf6<0)]=0;


plot(t_i,A_i,'go',label='Modeled Mean Age',ms=8);
plot(t_i,pist_age,'bo',label='3H/3He Piston Flow Age',ms=8);
plot(t_i,pist_age_cfc12,'rs',label='CFC12 Piston Flow Age',ms=8);
plot(t_i,pist_age_sf6,'k^',label='SF6 Piston Flow Age',ms=8);
ylabel('GW Age (yr)');
xlabel('year');
legend(loc='best',numpoints=1);
show();

