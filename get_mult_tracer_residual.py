"""
get_mult_tracer_residual.py - this script will print the objective function value for a 
dakota run.  Uses pflotran_tools.  
    Inputs - user needs to provide the directions to get d_obs
        the observed data.  This example is for a synthetic data set.  The correct
        file handle(s) for the observation points from the model output need to be input.
        This will be dependent upon the model, descritization, processor number.  Run the 
        initial model to get the observation point file names.  
    Outputs - Will print the value(s) of the objective function dependent upon the lst_sq flag.
        If lst_sq flag = 1 the difference at each point in time for each point one difference per line.
        If lst_sq_flag = 0 one line with the value of l2norm(Am-d).  This output should be directed 
        into the file which dakota will read in the simulator script for a dakota run.

Written W. Payton Gardner
"""


from pflotran_tools import *
from datetime import datetime

# get the observed data.
# for synthetic dataset 
#df_obs = ReadPflotranObsPtPandas(obsptfile='../1d_tracer_his.truth',date_start=datetime(1930,1,1));
df_obs = ReadPflotranObsPtPandas();  #for testing

# get the currently modeled data
#df_mod = ReadPflotranObsPtPandas(obsptfile='1d_tracer_his.dat',date_start=datetime(1930,1,1));
df_mod = ReadPflotranObsPtPandas();  #for testing

# now calculate the difference between the two
iVars=['P','He3_tot_M','He4_tot_M','CFC11_tot_M','CFC12_tot_M','CFC113_tot_M','SF6_tot_M','H3_tot_M'];

CalcResidualObsPt(df_obs,df_mod,iVars,lst_sq_flag=0);

