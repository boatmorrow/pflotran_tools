# -*- coding: utf-8 -*-
"""
pflotran tools.  A set of functions to help the python user deal with pflotran output 
Created on Fri Feb 24 15:27:45 2012

@author: W. Payton Gardner 
"""

from numpy import *
from scipy.integrate import trapz
import string
import pdb
from pylab import *
import os
import h5py
import pandas
import  datetime as dt
import matplotlib.dates as mdates

def mol_kg2tu(x):
    '''converts x molality to tu.'''
    molal2tu = 1/(1./10**18*(2.*6.02e23)/18./6.02e23*1000.);
    tu = float(x)*molal2tu;
    return tu

def mol_kg2ccSTP_g(x):
    '''convert x molality to ccSTP/g'''
    molal2cc = 22414.*1000;
    cc = float(x)*molal2cc;
    return cc

def read_pflotran_obs_pt_flo(ffile='observation-0.tec'):
    '''read_pflotran_obs_pt.py - Reads a pflotran observation point dataset (datafile) in tecplot point form and returns 
        a record array of the data DD.  For now only works for with flow data'''
        # I could do a more sophisticated/general naming process as below here, but for now this'll work
    dd = loadtxt(ffile,skiprows=1,dtype={'names':('time','Pres','Sl','vlx','vly','vlz'),'formats':('f4','f4','f4','f4','f4','f4')});
    return dd

def read_pflotran_obs_pt_conc(ffile='observation-0.tec'):
    '''read_pflotran_obs_pt.py - Reads a pflotran observation point dataset (datafile) in tecplot point form and returns 
        a record array of the data DD.  For now only works for with flow data'''
        # I could do a more sophisticated/general naming process as below here, but for now this'll work
    dd = loadtxt(ffile,skiprows=1,dtype={'names':('time','Conc'),'formats':('f4','f4')});
    return dd

def read_pflotran_output(tfile):
    '''read_pflotran_output(tfile) reads a pflotran tecplot point file and returns the time of the output, the data array dd'''
    df = open(tfile,'r');    
    lf = df.readlines();
    if len(string.split(lf[2]))>7:
        dt = float(string.split(lf[2])[3][0:-2]); #get the simulation time
    else:
        dt = float(string.split(lf[2])[1][3:-2]);
        
#get the variables output from the tecplot file
    p = string.split(lf[1],sep='"');
    p = array(p[1:-1]);
    p = p[where(p!=',')];
    format_list = [];
    name_list = list(p);
    #set_trace();
    for i in range(len(name_list)):
        format_list.insert(i,'f4');
        name_list[i] = string.strip(name_list[i]);
        try: 
            string.split(string.split(name_list[i])[0],'-')[1];
        except IndexError:
            name_list[i] = string.split(name_list[i])[0];
        else:
            name_list[i] = string.split(string.split(name_list[i])[0],'-')[1];
        
    format_list[-1] = 'i4';  #material id
    dd = loadtxt(tfile,skiprows=3, dtype = {'names':name_list, 'formats':format_list});
    return dt, dd

def plot_pflotran_output_line(dd,dt,direction,variable,pcolor='b',style='-'):
    ''' plots a pflotran variable over the direc`tion specified.  Needs the data array from reading in a pfltran output file
    and the sim time from the data file read.'''
    if direction == 'Z':
        plot(dd[variable],dd[direction],label='time = '+ str(dt),color=pcolor,linestyle=style);
        xlabel(variable);
        ylabel('distance ' + direction);
        #legend();
        #title('Time = ' + str(dt));
    else:
        plot(dd[direction],dd[variable],label='time = '+ str(dt),color=pcolor,linestyle=style);
        xlabel('distance ' + direction);
        ylabel(variable);
        #legend();
        #title('Time = ' + str(dt));
        
def plot_observation_pt_1d(dd):
    '''pretty simple/limited function to plot the concentration at a point.  only works for 1 tracer right now'''
    plot(dd['time'],dd['Conc']); #output from read_pflotran_obs_pt_conc
    xlabel('time(Yrs)'); #will be dependent upon problem
    ylabel('Concentration');

def get_evolution(dirpath='.'):
    '''makes a data dictionary keyed by time for all pflotran.tec files in the directory dirpath'''
    f_list = os.listdir(dirpath); # list of the files in dirpath
    d_dict = {};
    for f in f_list:
        ext = string.split(f,'.')[-1];
        prefix = string.split(string.split(f,'.')[0],sep='-')[0];
        if ext == 'tec':   # a tecplot file
            if prefix != 'observation':
                dt, dd = read_pflotran_output(f); # read the file
                d_dict[dt] = dd;  # add to dictionary
        else:
            continue
        
    return d_dict

def make_tracer_evo_plot(dd_dict):
    ''' Simple routine for making a tracer breakthrough plot'''
    for k in dd_dict.keys():
        dd = dd_dict[k];
        plot(dd['X'],dd['Tracer_tot_M'],label='time = '+ str(k));


def get_final_data(dd_dict):
    '''gets the final time and data dictionary'''
    max_t=amax(list(dd_dict.keys()));
    dde = dd_dict[max_t];
    return max_t,dde

def check_sum_max_spatial(dd):
    '''integrates the travel time pdf'''
    sum_prob = trapz(dd['Tracer_tot_M'],dd['X']);
    return sum_prob

def read_pflotran_h5(h5file='pflotran.h5'):
    '''returns arrays of x,y,z coordinates which give the corners, and data dictionary of the block centered values? 
    at each time'''
    f1 = h5py.File(h5file);
    x_n = f1['Coordinates']['X [m]'][:];
    y_n = f1['Coordinates']['Y [m]'][:];
    z_n = f1['Coordinates']['Z [m]'][:];
    #get block center coordinates
    
def ReadPflotranObsPtPandas(obsptfile='observation-59.tec',map_time2date=False,date_start=dt.datetime(1930,0o1,0o1)):
    '''reads in the observation point file from plotran and returns a pandas dataframe indexed by the 
    date.  Needs to know the start date of the model timesteps.  Assumes that the model time recorded
    by pflotran is date_start + time.'''
    ################################################################################
    # read in the output observation point
    ################################################################################
    f = file(obsptfile,'r');
    header = f.readline()
    #first we need to deal with that nast ass header
    lh = header.split(',')
    for i in range(len(lh)):
        lh_i = lh[i].strip();    
        lh_i = lh_i.strip('"');
        #lh[i]= string.join(lh_i.split(' ')[0:2],'_');  #just keep the first part of the original
        #lh[i]=lh_i.split(' ')[0]
    #now write the new mangle observation file (for now until I figure out to just change this to a Pandas data frame)
    print("reading in file")
    df = pandas.read_csv(f, sep='\s+', skiprows=0,names=lh,index_col=0)
    f.close()
    print("done reading, now mapping time")
    if map_time2date:
        try:
            days = df['Time_[y]']*365.25;  #change the plotran timestep output in years to days (this may change with different models - not sure...)
        except KeyError:
            days = df['Time_[yr]']*365.25
        years = ones(len(days));
        dyear = date_start;  
        df['Date']=dyear; #stat out at the start date
        
        for i in range(len(days)):
            d = dt.timedelta(days=days[i]);
            df['Date'][i] = df['Date'][i]+d;    # the true time is modeled time plus the start_time
        #set_trace();
        df_mod = df.set_index('Date');  # the pandas dataframe indexed on time
    else:
        df_mod=df
    return df_mod

def ReadPflotranObsPt(obsptfile='observation-59.tec'):
    '''reads in the observation point file from plotran and returns a numpy rec array.  
   With names given by the header in observation point file.'''
    ################################################################################
    # read in the output observation point
    ################################################################################    
    f = file(obsptfile,'r');
    ll = f.readlines();
    f.close();
    #first we need to deal with that nast ass header
    header = ll[0];
    lh = header.split(',');
    fm = [];
    ll2 = [];
    ll2.append(ll[0]);
    
    for line in ll[1:-1]:
        j=0;
        li = line.split(); 
        for field in li:
            if len(field.split('E'))==1:
                #set_trace();
                if len(field.split('-'))==2:
                    field = field.split('-')[0][0:-1] + 'E-' + field.split('-')[1];
                if len(field.split('+'))==2:
                    field = field.split('+')[0][0:-1] +'E+' + field.split('+')[1];    
                li[j]=field;
            j=j+1;
        line = string.join(li);
        line = line+'\n'
        ll2.append(line);
    manglefile = 'observation_pg.tec'
    f = file(manglefile,'w');
    f.writelines(ll2);
    f.close();
    #make the names and formats to define the rec array
    for i in range(len(lh)):
        lh_i = lh[i].strip();    
        lh_i = lh_i.strip('"');
        lh[i]= string.join(lh_i.split(' ')[0:2]);  #just keep the first part of the original
        #lh[i]=lh_i.split(' ')[0];
        fm.append('f4');
    
    dd = loadtxt(manglefile,skiprows=1,dtype={'names':lh,'formats':fm});
    return dd
    
def CalcResidualObsPt(df_obs,df_mod,iVars,lst_sq_flag=0):
    '''Calculates the residuals between observed and modeled output for a pflotran observation point.  df_mod
    and df_obs are modeled and observed data at the same point in the form of Pandas DataFrames keyed by the 
    time and date of sampling.  This script will find the closest previously modeled date to the observed date and compare
    and compare data for all variables in the list iVars.  iVars is a list of the column names which contain
    the data.  They must match for both modeled and observed data.  If lst_sq_flag is true then a the function will print
    the difference between the modeled and observed data at each point, one difference per line.  If lst_sq_flag
    is false then the l2 norm of the residual will be printed.'''
    ################################################################################
    # grab out the modeled data nearest the observation data and calculate the goal
    # function.  Should either be a vector of differences (for Dakota Least Squares),
    # or the l2 norm of the residual for the non-least squares methods.
    ################################################################################
    
    #here's one way to get the subset of modeled observations to match with the observations.  This will use the
    #previous modeled output in the case of a difference between model output times and observation times!
    df_i = df_mod.reindex(df_obs.index,columns=iVars,method='pad');
    #set_trace();
    #calcualte a least squares norm
    if lst_sq_flag:
        
        lines = array([]);
        for v in iVars:
            lines = concatenate((lines,df_i[v].values-df_obs[v].values));
            #may want to add a normalization routine here...
        
        for l in lines:
            print(l);  # print because in the dakota scriptology this output head to dakota.
    
    else:
        Am = array([]);
        d = array([]);
        for v in iVars:
            Am = concatenate((Am,df_i[v].values));
            d = concatenate((d,df_obs[v].values));
        print(linalg.norm((Am-d)));

def make_obs_pt_hist_plot(dd,iVar):
    '''Plots up a history of iVar vs. time for a given obs pt which has been read in to a Pandas data frame
    from ReadPflotranObsPtPandas. iVar is a list of variables'''
    #formatting for a date plot
    decades    = mdates.YearLocator(10);   # every year
    years   = mdates.YearLocator();  # every month
    yearsFmt = mdates.DateFormatter('%Y');
   
    #plot the data up
    fig = figure();
    ax = fig.add_subplot(111);
    for i in iVar:
        if i == 'H3_tot_M':
            ax.plot(dd.index,dd[i].map(mol_kg2tu),label=r'$^3$H (TU)');  
        if i == 'He3_tot_M':
            ax.plot(dd.index,dd[i].map(mol_kg2ccSTP_g),label=r'$^3$He (ccSTP/g)');  
        else:
            ax.plot(dd.index,dd[i],label=i);
    
    # format the ticks
    ax.xaxis.set_major_locator(decades);
    ax.xaxis.set_major_formatter(yearsFmt);
    ax.xaxis.set_minor_locator(years);
    
    #format the coords box
    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    ax.grid(True)
    xlabel('Year');
    ylabel('Modeled Value');
    legend();
    fig.autofmt_xdate()
    show();

def make_modern_recharge(df,model_time_0=dt.datetime(1776,1,1)):
    '''produced two files the list of constrains and the condition for modern recharge for each date (limited to daly right now), df needs to be a time indexed data frame in mol/kg e.g. produced by get_atm_conc, and the date where the model begins'''
    f = file('Modern_Constraints.txt','w')
    s = []
    for i in range(len(df)):
        s.append('CONSTRAINT ' + df.index[i].isodate() +'\n')
        s.append('  CONCENTRATIONS\n')
        for j in range(len(df.columns)):
            s.append('    '+df.columns[j]+'\t'+'%1.4E' %df.ix[i][j]+'\tF\n')
        #    s.append('    A(aq)'+'\t'+'%1.4E' %1.e-8+'\tF\n');  #for plotran water_age simulation
        s.append('    Tracer'+'\t'+'%1.4E' %1.e-16+'\tF\n');  #for plotran water_age simulation
        #    s.append('    Tracer_Age'+'\t'+'%1.4E' %1.e-16+'\tF\n');  #for plotran water_age simulation
    s.append('  /\n');
    s.append('END\n');
    s.append('\n')
    f.writelines(s);
    f.close();

    ff = file('Modern_Recharge_Condition.txt','w');
    s = [];
    s.append('TRANSPORT_CONDITION modern_recharge\n');
    s.append('  TYPE dirichlet_zero_gradient\n');
    s.append('  CONSTRAINT_LIST\n')
    for i in range(df):
        seconds = (df.index[i] - model_time_0).total_seconds #needs to be in seconds for pflotran
        sec_str = '%1.5E' %seconds
        const_nm = df.index[i].isotdate()
        s.append('    '+sec_str+'\t'+const_nm+'\n')
    s.append('  /\n');
    s.append('END');
    ff.writelines(s);
    ff.close();
