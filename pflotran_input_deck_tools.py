'''module containing utility for creating and manipulating pflotran input decks.'''
import pdb
import numpy as n
import string

class layer:
    '''the layer class should give the gridding all it needs.  the layer number is the zero indexed number of layer increasing w to e, thickness is the length of layer in given direction, the nn the number of nodes desired and the biased is a number code, 0 no bias, 1 single bias decreasing in the positive direction, 2 single bias increasing the positive direction and 3 dual bias.'''
    def __init__(self,number,thickness,nn,biased):
        self.num = number
        self.h = thickness
        self.n_step = nn
        self.bias = biased

class pflotran_strata_region:
    '''this class contains the data to define a region in a pflotran input deck. dim is the tuple (lenx,leny,lenz), loc is the tuple (x,y,z) of the lower right hand corner of the region, and strata_nm is the name of the strata to which the which the region should be coupled.'''
    def __init__(self,name,dim,loc,strata_nm,ic_cond_nm,cond_nm):
        self.name = name
        self.dim = dim
        self.loc = loc
        self.strata_nm = strata_nm
        self.ic_cond_nm = ic_cond_nm
        self.cond_nm = cond_nm

def write_region_card(reg,rlist):
    ''' writes a region card for a rectangular region class'''
    rlist.append('REGION ' + reg.name +'\n')
    rlist.append('  COORDINATES'+'\n');
    xmin = reg.loc[0];
    xmax = xmin+reg.dim[0];
    ymin = reg.loc[1];
    ymax = ymin+reg.dim[1];
    zmin = reg.loc[2];
    zmax = zmin+reg.dim[2];
    rlist.append('    '+str(xmin)+'  '+str(ymin)+'    '+str(zmin)+'\n');
    rlist.append('    '+str(xmax)+'  '+str(ymax)+'    '+str(zmax)+'\n');
    rlist.append('  /'+'\n');
    rlist.append('END'+'\n');
    rlist.append('\n');

def write_strata_card(reg,slist):
    '''writes a strata card for the region class'''
    # the strata cards for the stratigraphy
    slist.append('STRATA' +'\n');
    slist.append('  REGION ' + reg.name + '\n');
    slist.append('  MATERIAL ' + reg.strata_nm + '\n');
    slist.append('END\n');
    slist.append('\n');

def write_ic_condition_coupler_card(reg,clist,cname):
    '''write a condition copuler for a region class'''
    clist.append('INITIAL_CONDITION ' + cname +'\n');
    clist.append('  FLOW_CONDITION ' + reg.ic_cond_nm[0] +'\n');
    clist.append('  TRANSPORT_CONDITION ' + reg.ic_cond_nm[1] +'\n');
    clist.append('  REGION ' + reg.name +'\n');
    clist.append('END\n');
    clist.append('\n');

def write_source_sink_coupler_card(reg,sslist,ssname):
    '''write a source sink coupler card for a region class'''
    sslist.append('SOURCE_SINK ' + ssname +'\n');
    sslist.append('  FLOW_CONDITION ' + reg.cond_nm[0] +'\n');
    sslist.append('  TRANSPORT_CONDITION ' + reg.cond_nm[1] +'\n');
    sslist.append('  REGION ' + reg.name +'\n');
    sslist.append('END\n');
    sslist.append('\n');
def calc_grav_vec(dip,strike,g_mag=-9.8086):
    '''Calculate the coefficients of the gravity vector for a given slope angle and aspect -  g_hat = Ai_hat+Bj_hat+Ck_hat.
          Requires the slope angle in degrees and the strike angle where 0 is north. Returns the tuple (A,B,C).'''
    Fz = n.cos(dip*n.pi/180.)*g_mag #by definition of dip
    Fxy = n.cos((90.-dip)*n.pi/180.)*g_mag
    Fy = n.cos((90.-strike)*n.pi/180.)*Fxy
    Fx = n.cos(strike*n.pi/180.)*Fxy
    g = [Fx,Fy,Fz]
    for i in xrange(len(g)):
        if abs(g[i]) <= 1e-3:
            g[i] = 0.
    return tuple(g)

def write_grid_card_bounds(nxyz,b_min,b_max,grav=(0.,0.,0.),origin=(0.,0.,0.),invert_z=False):
    '''write a pflotran grid card for a rectanular grid.  
            
            Arguments:
                nxyz - the upple (nx, ny, nz)
                b_min - the tuple (xmin,ymin,zmin)
                b_max - the tuple (xmax,ymax,zmax)
                grav - the gravity vector (Fx,Fy,Fz)
                origin - coords of origin
                invert_z - z is positive downward.
    
            Returns:
                grid.txt - a file containing the grid card
       
      '''
    
    gc = []
    f = file('grid.txt','w')
    gc.append('GRID \n')
    gc.append('  TYPE structured \n')
    gc.append('  GRAVITY '+str(grav[0])+'  '+str(grav[1])+'  '+str(grav[2])+'\n')
    gc.append('  ORIGIN  '+str(origin[0])+'  '+str(origin[1])+'  '+str(origin[2])+'\n')
    if invert_z:
        gc.append('  INVERT_Z\n')
    gc.append('  NXYZ  '+str(int(nxyz[0]))+'  '+str(int(nxyz[1]))+'  '+str(int(nxyz[2]))+'\n')
    gc.append('  BOUNDS\n')
    gc.append('    '+str(b_min[0])+'  '+str(b_min[1])+'  '+str(b_min[2])+'\n')
    gc.append('    '+str(b_max[0])+'  '+str(b_max[1])+'  '+str(b_max[2])+'\n')
    gc.append('  /\n')
    gc.append('END')
    f.writelines(gc)
    print 'file written'
    return

def change_por(ifile,ofile,material,por):
    '''change the permeability in the material card.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'MATERIAL_PROPERTY'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == material
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                continue
            else:
                if line.split()[1] == material:
                    iflag = 1
                    continue

        if iflag:
            if line.split()[0][0:5] == 'POROS':
                ll[i] = '    '+ line.split()[0] +' ' + str(por) + '\n'
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

def change_perm(ifile,ofile,material,perm):
    '''change the permeability in the material card.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'MATERIAL_PROPERTY'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == material
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                continue
            else:
                if line.split()[1] == material:
                    iflag = 1
                    continue

        if iflag:
            if line.split()[0][0:5] == 'PERM_':
                ll[i] = '    '+ line.split()[0] +' ' + str(perm) + '\n'
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

def change_flow_condition_flux(ifile,ofile,bc_name,flux):
    '''change the the neuman flux.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    fflag = 0
    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'FLOW_CONDITION'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == bc_name
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                continue
            else:
                if line.split()[1] == bc_name:
                    iflag = 1
                    continue

        if iflag:
            if line.split()[0] == 'FLUX':
                if line.split()[1] =='neumann':
                    fflag=1
                    continue
        if fflag:
            if line.split()[0] =='FLUX':
                ll[i] = '  '+ line.split()[0] +' ' + str(flux) + '\n'
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

def change_flow_condition_rate(ifile,ofile,bc_name,rate):
    '''change the the volumetric injection rate.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    fflag = 0
    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'FLOW_CONDITION'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == bc_name
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                continue
            else:
                if line.split()[1] == bc_name:
                    iflag = 1
                    continue

        if iflag:
            if line.split()[0] == 'RATE':
                if line.split()[1] =='volumetric_rate':
                    fflag=1
                    continue
        if fflag:
            if line.split()[0] =='RATE':
                ll[i] = '  '+ line.split()[0] +' ' + str(rate) + '\n'
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

def change_cap_press_params_vg(ifile,ofile,sat_func,m,alpha):
    '''change the capillary pressure parameters for the van genuchten model.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    lflag = 0

    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'CHARACTERISTIC_CURVES'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == sat_func
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                    lflag = 0
                continue
            else:
                if line.split()[0] == 'CHARACTERISTIC_CURVES':
                    if line.split()[1] == sat_func:
                        lflag = 1
                        continue
        
        if lflag:
            if line.split()[0] == 'SATURATION_FUNCTION':
                iflag = 1
                continue
            if line.split()[0] == 'END':
                lflag = 0
                continue

        if iflag:
            if line.split()[0] == 'M':
                ll[i] = '    '+ line.split()[0] +' ' + str(m) + '\n'
            if line.split()[0] == 'ALPHA':
                ll[i] = '    '+ line.split()[0] +' ' + str(alpha) + '\n'
            if line.split()[0] == 'PERMEABILITY_FUNCTION':
                iflag = 0
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

def change_rel_perm_params_mualem(ifile,ofile,sat_func,m,res_sat):
    '''change the relative permeability parameters for the van genuchten-mualem model.'''

    #open the inputfile
    f = open(ifile,'r')
    ll=f.readlines()
    f.close()

    #parse and replace relavent fields
    iflag = 0
    lflag = 0

    for i in xrange(len(ll)):
        line = ll[i]
        try:
            line.split()[0] == 'CHARACTERISTIC_CURVES'
        except IndexError:
            continue
        else:
            try:
                line.split()[1] == sat_func
            except IndexError:
                if line.split()[0] == 'END':
                    iflag = 0
                    lflag = 0
                continue
            else:
                if line.split()[0] == 'CHARACTERISTIC_CURVES':
                    if line.split()[1] == sat_func:
                        lflag = 1
                        continue
        
        if lflag:
            if line.split()[0] == 'PERMEABILITY_FUNCTION':
                iflag = 1
                continue
            if line.split()[0] == 'END':
                lflag = 0
                continue

        if iflag:
            if line.split()[0] == 'M':
                ll[i] = '    '+ line.split()[0] +' ' + str(m) + '\n'
            if line.split()[0] == 'LIQUID_RESIDUAL_SATURATION':
                ll[i] = '    '+ line.split()[0] +' ' + str(res_sat) + '\n'
            if line.split()[0] == 'END':
                iflag = 0
    
    #write the new file
    f = open(ofile,'w')
    f.writelines(ll)
    f.close()
    return

#should do it for a list of patterns and cardfiles...
def insert_card(ifile,ofile,pattern,cardfile):
    '''inserts the card in cardfile at location pattern in ifile.  returns ofile'''
    f = file(ifile,'r')
    ll = f.readlines()
    f.close()

    f = file(cardfile,'r')
    il = f.readlines()
    f.close()
    
    #now march through and find the right places to insert
    for i in xrange(len(ll)):
        line = string.strip(ll[i])
        if line == pattern:
            for x in reversed(il):
                ll.insert(i+1,x)  
    
    f = file(ofile,'w')
    f.writelines(ll)
    f.close()
