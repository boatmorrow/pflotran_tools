'''module containing utility for creating and manipulating pflotran input decks.'''
import pdb
import numpy as n
import string
from scipy import optimize


class layer:
    '''the layer class should give the gridding all it needs.  the layer number is the zero indexed number of layer increasing w to e, thickness is the length of layer in given direction, the nn the number of nodes desired and the biased is a number code, 0 no bias, 1 single bias decreasing in the positive direction, 2 single bias increasing the positive direction and 3 dual bias.'''
    def __init__(self,number,thickness,nn,biased,name='Who Dey'):
        self.num = number
        self.h = thickness
        self.n_step = nn
        self.bias = biased
        self.name = name

class pflotran_strata_region:
    '''this class contains the data to define a region in a pflotran input deck. dim is the tuple (lenx,leny,lenz), loc is the tuple (x,y,z) of the lower right hand corner of the region, and strata_nm is the name of the strata to which the which the region should be coupled.'''
    def __init__(self,name,dim,loc,strata_nm,ic_cond_nm,cond_nm):
        self.name = name
        self.dim = dim
        self.loc = loc
        self.strata_nm = strata_nm
        self.ic_cond_nm = ic_cond_nm
        self.cond_nm = cond_nm


class pflotran_bc_region:
    '''this class contains the data to define a boundary face region in a pflotran input deck. dim is the tuple (lenx,leny,lenz), one of which should be zero, loc is the tuple (x,y,z) of the lower right hand corner of the region, and face_nm is the name of the face north, south, east, west, top, bottom.'''
    def __init__(self,name,dim,loc,face_nm,cond_nm):
        self.name = name
        self.dim = dim
        self.loc = loc
        self.face_nm = face_nm
        self.cond_nm = cond_nm

class pflotran_observation_point:
    '''this class contains information for a pflotran observation point, and the region associated with it.'''
    def __init__(self,name,loc):
        self.name = name
        self.loc = loc
    def write_observation_card(self,olist):
        '''write the observation card'''
        olist.append('OBSERVATION'+'\n')
        olist.append('  REGION '+self.name+'\n')
        olist.append('  VELOCITY'+'\n')
        olist.append('END'+'\n')
        olist.append('\n')
        return olist
    def write_region_card(self,rlist):
        '''write the region card'''
        rlist.append('REGION ' + self.name + '\n')
        rlist.append('  COORDINATE '+str(self.loc[0])+' '+str(self.loc[1])+' '+str(self.loc[2])+'\n')
        rlist.append('END'+'\n')
        rlist.append('\n')
        return rlist

def write_region_card(reg,rlist):
    ''' writes a region card for a rectangular region class. returns rlist and list of the strings with the region card. '''
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
    return rlist

def write_region_card_bc(reg,rlist):
    ''' writes a region card for a face/bc region class. returns rlist and list of the strings with the region card. '''
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
    rlist.append('  FACE '+reg.face_nm+'\n')
    rlist.append('END'+'\n');
    rlist.append('\n');
    return rlist
    
def write_strata_card(reg,slist):
    '''writes a strata card for the region class'''
    # the strata cards for the stratigraphy
    slist.append('STRATA' +'\n');
    slist.append('  REGION ' + reg.name + '\n');
    slist.append('  MATERIAL ' + reg.strata_nm + '\n');
    slist.append('END\n');
    slist.append('\n');

def write_ic_condition_coupler_card(reg,clist,cname,tran_flag=True):
    '''write a condition copuler for a region class'''
    clist.append('INITIAL_CONDITION ' + cname +'\n');
    clist.append('  FLOW_CONDITION ' + reg.ic_cond_nm[0] +'\n');
    if tran_flag:
        clist.append('  TRANSPORT_CONDITION ' + reg.ic_cond_nm[1] +'\n');
    clist.append('  REGION ' + reg.name +'\n');
    clist.append('END\n');
    clist.append('\n');

def write_source_sink_coupler_card(reg,sslist,ssname,tran_flag=True):
    '''write a source sink coupler card for a region class'''
    sslist.append('SOURCE_SINK ' + ssname +'\n');
    sslist.append('  FLOW_CONDITION ' + reg.cond_nm[0] +'\n')
    if tran_flag:
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
    f.close()
    print 'file written'
    return

def write_layered_dx(x):
    '''dxx = write_layered_dx(x)
        Returns dxx - a list of grid distance increments for a series of layers with set thicknesses, desired nodes, and bias flag (see layer class).
        Requires x - a list of layer classes with defined thickness number of nodes and bias type.  Layers should be ordered in pflotran natural ordering - (i.e. starting a lower southwest corner.'''
    dx = n.array([])
    l_count = 0
    txx = 0.
    for l in x:
        txx = txx + l.h
        bias_l = l.bias
        nx = l.n_step
        x1 = l.h/l.n_step # just a guess for now
        if l.bias == 0:
            px = sum(dx)
            tx = px + l.h
            for i in xrange(nx):
                dx = n.append(dx,x1)
            #hit end of layer on the nose
            dx[-1] = tx - sum(dx[0:-1])
        if l.bias == 1:  #single bias decreasing in positive direction
            if x[l_count+1].bias == 0:
                x1 = x[l_count+1].h/x[l_count+1].n_step # want to end with the next grid spacing
            else:
                print 'not sure what do with bias min x distance, figure it out'
            args=(l.h,nx,x1)
            bx = optimize.bisect(txb_diff_single_bias,1.,2.5,args=args);
            px = sum(dx)
            tx = px+l.h
            #construct the biased vector
            for i in xrange(1,int(nx)):
                dx = n.append(dx,x1*bx**(nx-(i+1)))
            #hit end of bias on the nose
            dx[-1] = tx - sum(dx[0:-1])
        if l.bias == 2: #single bias increasing in positive direction
            if l.num >= 1:
                x1 = dx[-1] #start with last grid spacing
            else:
                print 'not sure what to do about min x distance, figure it out'
            bx = optimize.bisect(txb_diff_single_bias,1.,2.5,args=(l.h,nx,x1));
            px = sum(dx)
            tx = px+l.h
            #construct the biased vector
            for i in xrange(1,int(nx)):
                dx = n.append(dx,x1*bx**(i-1));

            #hit end of domain on the nose
            dx[-1] = tx - sum(dx[0:-1])
        if l.bias == 3: #dual bias
            if l.num >= 1:
                x1 = dx[-1] #start with last grid spacing
            else:
                print 'not sure what to do about min x distance, figure it out'
            # now to figure out the correct bias (for a dual bias)
            nnx = nx/2-1;
            px = sum(dx)
            tx = px+l.h
            bx = optimize.bisect(txb_diff_dual_bias,1.,2.5,args=(l.h,nx,x1));
            #first half of the biased vector
            for i in xrange(1,int(nx/2)):
                dx = n.append(dx,x1*bx**(i-1));
            #second half of the biased x vector
            for i in xrange(1,int(nx/2)):
                dx = nappend(dx,x1*bx**(nnx-(i)));
            #hit end of bias on the nose
            dx[-1] = tx - sum(dx[0:-1]);
        tl = n.sum(dx)-px
        print 'Total thickness of layer ' + l.name + ' = %1.2g units' % tl
        l_count += 1

    #make the ending length get us to a the desired total length.
    print "layer discretization error %1.5g" % (txx-tx)
    dx[-1] = txx - sum(dx[0:-1])

    dxx=list(dx)
    print 'Total layered thickness' + l.name +' = %1.2g units' %(sum(dxx))
    return dxx

def write_grid_card_dx(dx,dy,dz,grav=(0.,0.,0.),origin=(0.,0.,0,),invert_z=False):
    '''write a pflotran grid card for a rectanular grid using incremental distances.
            
            Arguments:
                dx,dy,dz = array like increments in the x,y,z direction
                grav - the gravity vector (Fx,Fy,Fz)
                origin - coords of origin
                invert_z - z is positive downward.

            Returns:
                grid.txt - a file containing the grid card
       
      '''
    #print out the grid card
    dx_len=str(len(dx));
    dz_len=str(len(dz));
    dy_len=str(len(dy));
    gfile = file('grid.txt','w');
    i=0;
    for x in dx:
        if i%10.==0:
            if i != 0:
                dx.insert(i,'\\'+'\n')
                i=i+1;
                continue
        dx[i]='%13.11e' %dx[i]+' '
        i=i+1

    i=0
    for y in dy:
        if i%10.==0:
            if i != 0:
                dy.insert(i,'\\'+'\n')
                i=i+1
                continue
        dy[i]='%13.11e' %dy[i]+' '
        i=i+1

    i=0;
    for z in dz:
        if i%10.==0:
            if i != 0:
                dz.insert(i,'\\'+'\n')
                i=i+1;
                continue
        dz[i]='%13.11e' %dz[i]+' '
        i=i+1

    dx.append('\n')
    dy.append('\n')
    dz.append('\n')

    gfile.write('GRID\n')
    gfile.write('  TYPE structured\n')
    gfile.write('  NXYZ '+dx_len+' '+dy_len+' '+dz_len+'\n')
    gfile.write('  DXYZ\n')

    gfile.writelines(dx)
    gfile.writelines(dy)
    gfile.writelines(dz)
    gfile.write('  /\n')
    gfile.write('END\n')
    gfile.close()
    print 'file written'

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

def txb_diff_single_bias(bx,tfx,nx,x1):
    '''returns the difference between the desired distance (args[0]) and the calculated distance for a bias of bx, with args[1] steps.'''
    tbx = 0;
    for i in xrange(1,(int(nx))):
        tbx = tbx+(x1*bx**(i-1));
    return tbx-tfx;


def txb_diff_dual_bias(bx,tfx,nx,x1):
    '''returns the difference between the desired distance (args[0]) and the calculated distance for a bias of bx, with args[1] steps. Dual biased vector'''
    tbx = 0;
    tfx = tfx/2;
    for i in xrange(1,(int(nx/2))):
        tbx = tbx+(x1*bx**(i-1));
    return tbx-tfx;
