# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 12:00:15 2012

@author: wpgardn
"""

import h5py
import pdb
from pylab import *
from numpy import *
import cfc_tools as cfc
import sf6_tools as sf6
import noble_gas_tools as ng
from matplotlib import ticker, cm, colors


# get the coordinate array turn them to cell center coordinates - I think this will only work for block grids...
def GetCellCenters(inputfile):
    '''Gets the cell centers for rectangular grids from a pflotran h5 file.  Returns xx,yy,zz the cell center
    arrays.'''
    f1 = h5py.File(inputfile,'r');
    #x
    x_n = f1['Coordinates']['X [m]'][:];
    xmin = x_n[0:len(x_n)-1];
    xmax = x_n[1:len(x_n)];
    xx = (xmin+xmax)/2;
    #y
    y_n = f1['Coordinates']['Y [m]'][:];
    ymin = y_n[0:len(y_n)-1];
    ymax = y_n[1:len(y_n)];
    yy = (ymin+ymax)/2;
    #z
    z_n = f1['Coordinates']['Z [m]'][:];
    zmin = z_n[0:len(z_n)-1];
    zmax = z_n[1:len(z_n)];
    zz = (zmin+zmax)/2;
    f1.close();
    return xx,yy,zz

#get time information
def GetTimeInfo(inputfile):
    '''Gets the time information and time keys for a pflotran output file.  Returns the output time in years, and the timekeys
    which is a list of of all the keys for different out put times in a given h5 file.'''
    f1 = h5py.File(inputfile,'r');    
    times = [];
    timekeys = [];
    #tracer = 'H3'
    #pdb.set_trace();
    for k in f1.keys():
        kk = k.split()[0];
        if kk =='Time:':
            times.append(float(k.split()[1]));
            timekeys.append(k);
            #he4_f = reshape(f1[k][tracer+'_tot_M'],(len(xmid)));
            #he4_f = reshape(f1[k]['CFC11_tot_M'],(len(xmid)));
            #plot(xmid,he4_f,label='t = '+k.split()[1]);
    #pdb.set_trace();
    #timekeys = sort(timekeys);
    #times = sort(times);
    timekeys = array(timekeys);
    times = array(times);
    f1.close();
    return times,timekeys

#x is column and y is row in the plotran h5

def GetMaxTime(inputfile):
    '''return the maximum time output from the input file'''
    d = hdf5.File(inputfile,'r')
    kl = max([float(k.split()[1]) for k in d.keys() if k.split()[0]==u'Time:'])
    return kl

def PlotHeadDistribution2dXY(time,z,inputfile):
    '''Plots the pressure head distribution for a 2-d isothermal, fresh-water slice in the XY plane at location z.  Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    zindex = (abs(zz-z)).argmin();
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    hxy = f1[tk]['Liquid_Pressure [Pa]'][:,:,zindex]/1000./9.8 + zz[zindex]
    Vmin = hxy.min()
    Vmax = hxy.max()
    V = linspace(Vmin,Vmax,50); #linspace for now.  could add logspace option...
    cs = contourf(xx,yy,transpose(hxy),V); #The flip flop is because of the x=column y=row
    tticks = linspace(V.min(),V.max(),10);
    bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label('Head');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    title('Time = '+str(time));
    #show();
    f1.close();
     
def PlotHeadDistribution2dXZ(time,y,inputfile):
    '''plots the pressure head distribution for a 2-d isothermal, fresh-water slice in the XZ plane at location y.  will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    yindex = (abs(yy-y)).argmin();
    f1 = h5py.File(inputfile,'r');
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    hxz = f1[tk]['Liquid_Pressure [Pa]'][:,yindex,:]/1000./9.8 + zz
    Vmin = hxz.min()
    Vmax = hxz.max()
    V = linspace(Vmin,Vmax,50); #linspace for now.  could add logspace option...
    cs = contourf(xx,zz,transpose(hxz),V); #The flip flop is because of the x=column y=row
    tticks = linspace(V.min(),V.max(),10);
    bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label('Head');
    xlabel('X distance (m)');
    ylabel('Z distance (m)');
    title('Time = '+str(time));
    #show();
    f1.close();
     
def PlotHeadDistribution2dYZ(time,x,inputfile):
    '''Plots the pressure head distribution for a 2-d isothermal, fresh-water slice in the yz plane at location x.  Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    xindex = (abs(xx-x)).argmin()
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    hyz = f1[tk]['Liquid_Pressure [Pa]'][xindex,:,:]/1000./9.8 + zz
    Vmin = hyz.min()
    Vmax = hyz.max()
    V = linspace(Vmin,Vmax,50); #linspace for now.  could add logspace option...
    cs = contourf(yy,zz,transpose(hyz),V); #The flip flop is because of the x=column y=row
    tticks = linspace(V.min(),V.max(),10);
    bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label('Head');
    xlabel('Y distance (m)');
    ylabel('Z distance (m)');
    title('Time = '+str(time));
    #show();
    f1.close();
     
def PlotMaterialDistribution2dXZ(time,y,inputfile):
    '''Plots the pressure material distribution for a 2-d xz slice at location y (i.e. x-sec).  Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile)
    #this is in model time 
    f1 = h5py.File(inputfile,'r')
    yindex = (abs(yy-y)).argmin()
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    mxz = f1[tk]['Material_ID'][:,yindex,:]
    Vmin = mxz.min()
    Vmax = mxz.max()
    V = arange(Vmin-1,Vmax+1,1.); #linspace for now.  could add logspace option...
    cs = contourf(xx,zz,transpose(mxz),V); #The flip flop is because of the x=column y=row
    #tticks = linspace(V.min(),V.max(),10);
    bar = colorbar()
    #bar.set_label('Head');
    xlabel('Y distance (m)');
    ylabel('Z distance (m)');
    title('Material Distribution');
    #show();
    f1.close();
     
def PlotVelocity2d(time,inputfile):
    '''Plots the velocity vectors for each cell in the bottom most row. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    #x component of vector
    U = transpose(f1[tk]['Liquid X-Velocity'][:,:,0]); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    #x component of vector
    V = transpose(f1[tk]['Liquid Y-Velocity'][:,:,0]); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    imshow(transpose(f1[tk]['Liquid_Pressure [Pa]'][:,:,0]/1000./9.8),origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label('Head (m)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    quiver(xx,yy,U,V);
    #show();
    f1.close();

#plot distriubtion of a tracer for a given time in xy plane
def PlotTracerDistribution2dXY(tracer,time,z,inputfile,logflag=False):
    '''Plots the tracer distribution in the xy plane for a given time at location z.  Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    tracerkey = 'Total_'+tracer+' [M]'
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    zi = abs(z-zz).argmin()
    tk = timekeys[(abs(time-times)).argmin()];
    Vmin = f1[tk][tracerkey][:,:,zi].min()
    Vmax = f1[tk][tracerkey][:,:,zi].max()
    V = linspace(Vmin,Vmax/10.,50); #linspace for now.  could add logspace option...
    if logflag:
        V_exp = arange(floor(log10(Vmin)-1),ceil(log10(Vmax)+1))
        V = power(10,V_exp)
        cs = contourf(xx,yy,transpose(f1[tk][tracerkey][:,:,zi]),V,norm=colors.LogNorm()); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    else:
        cs = contourf(xx,yy,transpose(f1[tk][tracerkey][:,:,zi]),V); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    #imshow(transpose(f1[tk][tracer][:,:,0]),origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
#    bar = colorbar(orientation='horizontal');
    tticks = linspace(V.min(),V.max(),10);
    if logflag:
        bar = colorbar(cs)
    else:
        bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label(tracer);
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    title('Time = '+str(time));
    #show();
    f1.close();

#plot distriubtion of a tracer for a given time in the xz plane
def PlotTracerDistribution2dXZ(tracer,time,y,inputfile,logflag=False):
    '''Plots the tracer distribution for a given time in the xz plane at location y. Will use the time key closest to the desired time.'''
    sxz,xx,zz = GetTracerXZSlice(tracer,y,time,inputfile)
    Vmin = sxz.min()
    Vmax = sxz.max()
    V = linspace(Vmin,Vmax/10.,50); #linspace for now.  could add logspace option...
    if logflag:
        V_exp = arange(floor(log10(Vmin)-1),ceil(log10(Vmax)+1))
        V = power(10,V_exp)
        cs = contourf(xx,zz,transpose(sxz),V,norm=colors.LogNorm()); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    else:
        cs = contourf(xx,zz,transpose(sxz),V); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    #imshow(transpose(f1[tk][tracer][:,:,0]),origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
#    bar = colorbar(orientation='horizontal');
    tticks = linspace(V.min(),V.max(),10);
    if logflag:
        bar = colorbar(cs)
    else:
        bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label(tracer);
    xlabel('X distance (m)');
    ylabel('Z distance (m)');
    title('Time = '+str(time));
    #show();

#plot distriubtion of a tracer for a given time in the yz plane
def PlotTracerDistribution2dYZ(tracer,time,x,inputfile,logflag=False):
    '''Plots the tracer distribution for a given time in the yz plane at location x. Will use the time key closest to the desired time.'''
    syz,yy,zz = GetTracerYZSlice(tracer,x,time,inputfile)
    Vmin = syz.min()
    Vmax = syz.max()
    V = linspace(Vmin,Vmax/10.,50); #linspace for now.  could add logspace option...
    if logflag:
        V_exp = arange(floor(log10(Vmin)-1),ceil(log10(Vmax)+1))
        V = power(10,V_exp)
        cs = contourf(yy,zz,transpose(syz),V,norm=colors.LogNorm()); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    else:
        cs = contourf(yy,zz,transpose(syz),V); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    #imshow(transpose(f1[tk][tracer][:,:,0]),origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
#    bar = colorbar(orientation='horizontal');
    tticks = linspace(V.min(),V.max(),10);
    if logflag:
        bar = colorbar(cs)
    else:
        bar = colorbar(ticks=tticks,format='%4.1e');
    bar.set_label(tracer);
    xlabel('Y distance (m)');
    ylabel('Z distance (m)');
    title('Time = '+str(time));
    #show();

def PlotMeanAgeDistribution(time,inputfile):
    '''Plots the tracer distribution for a given time. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    imshow(transpose(f1[tk]['Total Tracer_Age [M]'][:,:,0]/3.156e7),origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label('Mean Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def Plot3He3HAgeDistribution(time,inputfile):
    '''Plots the tritium helium age distribution for a given time. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    lambda_h3=0.0558;
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    c_he = transpose(f1[tk]['Total He3 [M]'][:,:,0]);
    c_h = transpose(f1[tk]['Total H3 [M]'][:,:,0]);
    # for now hard coded - but should be done smater in the future
    c_he_atm = 2.8253e-15; #taken from pflotran input file - sea level 25 C - zero salinty
    pist_age = 1./lambda_h3*log((c_he-c_he_atm)/c_h+1.);
    imshow(pist_age,origin='lower',extent=(0,max(xx),0,max(yy)),vmin=0,vmax=75); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label(r'$^3$H/$^3$He Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def PlotCFC12AgeDistribution(time,inputfile,t_start=1930):
    '''Plots the CFC12 age distribution for a given time. Assumes a hard coded water temperature of 25. C and elevation of zero m. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    t = times[(abs(time-times)).argmin()]+t_start; #move to date year
    #figure();
    c_cfc = transpose(f1[tk]['Total CFC12 [M]'][:,:,0]);
    #pdb.set_trace()
    vfunc = vectorize(cfc.equil_air_conc);
    vfunc2 = vectorize(cfc.age_date);
    eac = vfunc(c_cfc,25.0,0.0)/1e-12; #pptv
    pist_year_cfc12 = vfunc2(eac);
    pist_age_cfc12 = t-pist_year_cfc12;
    imshow(pist_age_cfc12,origin='lower',extent=(0,max(xx),0,max(yy)),vmin=0,vmax=75); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label('CFC-12 Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def PlotSF6AgeDistribution(time,inputfile,t_start=1930):
    '''Plots the SF6 age distribution for a given time. Assumes a hard coded water temperature of 25C and elevation of zero. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    t = times[(abs(time-times)).argmin()]+t_start; #move to date year
    #figure();
    c_sf6 = transpose(f1[tk]['Total SF6 [M]'][:,:,0]);
    #pdb.set_trace()
    vfunc = vectorize(sf6.equil_air_conc);
    vfunc2 = vectorize(sf6.age_date);
    eac = vfunc(c_sf6,25.0,0.0)/1e-15; #pptv
    pist_year_sf6 = vfunc2(eac);
    pist_age_sf6 = t-pist_year_sf6;
    imshow(pist_age_sf6,origin='lower',extent=(0,max(xx),0,max(yy)),vmin=0,vmax=75); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label(r'SF$_6$ Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def Plot39ArAgeDistribution(time,inputfile):
    '''Plots the 39Ar age distribution for a given time. Assumes a hard coded intial Ar39 concentration. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    c_Ar39 = transpose(f1[tk]['Total Ar39 [M]'][:,:,0]);
    lambda_Ar39=2.58e-3;
    c_Ar39_atm = 1.0211e-20;  # hard coded for now, but should be a function of recharge temp etc...
    pist_age_Ar39 = -lambda_Ar39**-1*log(c_Ar39/c_Ar39_atm);
    #pdb.set_trace();
    pist_age_Ar39[where(pist_age_Ar39<100)]=0;
    pist_age_Ar39[where(pist_age_Ar39>1000)]=1000;
    #Ar39_age.append(mean(pist_age_Ar39));
    #Ar39_std.append(std(pist_age_Ar39));
    #print 'Average 3H/3He Piston Flow Age = '+str(mean(pist_age)) + ' years'
    #print 'Standard Deviation of 3H/3He Piston Flow Age is '+str(std(pist_age))+ ' years'
    imshow(pist_age_Ar39,origin='lower',extent=(0,max(xx),0,max(yy))); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label(r'$^{39}$Ar Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def Plot81KrAgeDistribution(time,inputfile):
    '''Plots the 81Kr age distribution for a given time. Assumes a hard coded initial Kr81 concentration. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    c_Kr81 = transpose(f1[tk]['Total Kr81 [M]'][:,:,0]);
    lambda_Kr81=3.03e-6;
    c_Kr81_atm = 2.4938E-22;  # hard coded for now, but should be a function of recharge temp etc...
    pist_age_Kr81 = -lambda_Kr81**-1*log(c_Kr81/c_Kr81_atm);
    #pdb.set_trace();
    pist_age_Kr81[where(pist_age_Kr81<7.e4)]=0;
    pist_age_Kr81[where(pist_age_Kr81>7.e5)]=7.e5;   
    imshow(pist_age_Kr81,origin='lower',extent=(0,max(xx),0,max(yy)),vmin=0.,vmax=500.); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label(r'$^{81}$Kr Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();

def Plot4HeAgeDistribution(time,inputfile):
    '''Plots the 4He age distribution for a given time. Assumes a hard coded average crustal heliume production rate, recharge temp. of 25. and elevation of zero. Depth value is hard coded to the bottom layer right now. Will use the time key closest to the desired time.'''
    xx,yy,zz = GetCellCenters(inputfile);
    #this is in model time 
    f1 = h5py.File(inputfile,'r');
    times,timekeys = GetTimeInfo(inputfile);
    tk = timekeys[(abs(time-times)).argmin()];
    #figure();
    c_He4 = transpose(f1[tk]['Total He4 [M]'][:,:,0]);
    He4_prod = 5.756e-21;  #hard wired from the pflotran input file (mol/kg_h2o/sec - at some point should take advantage of my He tools suite.)
    c_He4_atm = ng.equil_conc('He',25.);  # hard coded for now, but should be a function of recharge temp etc...
    pist_age_He4 = (c_He4-c_He4_atm)/He4_prod;
    #pdb.set_trace();
    pist_age_He4[where(pist_age_He4<7.e4)]=0;
 
    imshow(pist_age_He4,origin='lower',extent=(0,max(xx),0,max(yy)),vmin=0.,vmax=500.); #z slice hardcoded for now, the flip flop is because of the x=column y=row
    bar = colorbar(orientation='horizontal');
    bar.set_label(r'$^{4}$He Piston Flow Age (yr)');
    xlabel('X distance (m)');
    ylabel('Y distance (m)');
    #show();
    f1.close();


#plot concentration at a point through time (i.e. obs pt)
def GetTracerObsSeries(tracer,x,y,z,inputfile):
    '''Plot the tracer concentration at a point for all times in the snapshot file.'''
    xx,yy,zz = GetCellCenters(inputfile);
    f1 = h5py.File(inputfile,'r');
    xindex = (abs(xx-x)).argmin();
    yindex = (abs(yy-y)).argmin();
    zindex = (abs(zz-z)).argmin();
    times,timekeys = GetTimeInfo(inputfile);
    ts = [];
    for k in timekeys:
        ck = f1[k][tracer][xindex,yindex,zindex];  #not sure if z index is correct yet
        #pdb.set_trace();
        ts.append(ck);
    #this is in model time 
    f1.close();
    return array(times),array(ts)

#get a y slice
def GetTracerXZSlice(tracer,y,time,inputfile):
    '''Returns an array of the tracer concentration in an xz plane at the y location at the time 
    t (in model time for now...)'''
    xx,yy,zz = GetCellCenters(inputfile);
    times,timekeys = GetTimeInfo(inputfile);
    f1 = h5py.File(inputfile,'r');
    yindex = (abs(yy-y)).argmin();
    tk = timekeys[abs(time-array(times)).argmin()];
    tracerkey = 'Total_'+tracer+' [M]'
    sxz = f1[tk][tracerkey][:,yindex,:];
    f1.close();
    return sxz,xx,zz

#get a x profile
def GetTracerXProfile(tracer,y,z,time,inputfile):
    '''Returns an array of the tracer concentration along an x progile at the y and z location at the nearest timekey
    in the model output file.'''
    tracerkey = 'Total_'+tracer+' [M]'
    xx,yy,zz = GetCellCenters(inputfile)
    times,timekeys = GetTimeInfo(inputfile)
    f1 = h5py.File(inputfile,'r')
    yindex = (abs(yy-y)).argmin()
    zindex = (abs(zz-z)).argmin()
    tk = timekeys[abs(time-array(times)).argmin()]
    sx = f1[tk][tracerkey][:,yindex,zindex]
    f1.close()
    return sx,xx
    
#get a x profile
def GetHeadXProfile(y,z,time,inputfile):
    '''Returns an array of the head along an x profile at the y and z location at the nearest timekey
    in the model output file.'''
    xx,yy,zz = GetCellCenters(inputfile)
    times,timekeys = GetTimeInfo(inputfile)
    f1 = h5py.File(inputfile,'r')
    yindex = (abs(yy-y)).argmin()
    zindex = (abs(zz-z)).argmin()
    tk = timekeys[abs(time-array(times)).argmin()]
    sx = f1[tk]['Liquid_Pressure [Pa]'][:,yindex,zindex]/1000./9.8+zz[zindex]
    f1.close()
    return sx,xx
    
#get a y slice
def GetTracerYZSlice(tracer,x,time,inputfile):
    '''Returns an array of the tracer concentration in an xz plane at the y location at the time 
    t (in model time for now...)'''
    xx,yy,zz = GetCellCenters(inputfile);
    #pdb.set_trace();
    times,timekeys = GetTimeInfo(inputfile);
    f1 = h5py.File(inputfile,'r');
    xindex = (abs(xx-x)).argmin();
    tk = timekeys[abs(time-array(times)).argmin()];
    tracerkey = 'Total_'+tracer+' [M]'
    sxz = f1[tk][tracerkey][xindex,:,:];
    f1.close();
    return sxz,yy,zz


#get the first and second moments of the mean age at many slices along the flow direction
def get_age_dist_info_w_dist(ifile,x_dist_vec,yr_start,output_yr):
    '''return the first and second moments of the age distribution, and the tracer derived ages at a y-z slice for each x
    in the x_distance vector.'''
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
    for x_dist in x_dist_vec:
        sxz,yy,zz = GetTracerYZSlice('Total Tracer_Age [M]',x_dist,output_yr-yr_start,ifile);   
        sxz = sxz/3.156e7;
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
    
    return MM_age, MM_std, H3He3M_age, H3He3_std, CFC12M_age, CFC12_std, SF6_age, SF6_std, Ar39_age, Ar39_std
