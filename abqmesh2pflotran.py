import h5py
import numpy
import pdb

inp_file = 'int_drift_mesh.inp';
h5_filename = 'mesh_int_drift.h5'
h5file = h5py.File(h5_filename,mode='w');
#element type number 8 is hex
elem_num = 8;

# open and parse the inp file

header = 0;

ifile = open(inp_file);
ll = ifile.readlines();
ifile.close();

# The way Cubit exports a sideset to Abaqus is with the entire element.
# So for plotran we need the sideset and the nodeset for a particular suface
# in order to define the faces like pflotran wants.  I'll make a dictionary of
# all vertices, and then a dictionary of all elements as I read the file.
# For very large meshes this will be a problem...

v_dict={};
e_dict={};
ss_count = 0;
ns_count = 0;
ns_flag = 0;

#vertex array
print "processing .inp file"
for line in ll:
   
    if line == "":
        continue;

    if line[0:4] == "**\n":
#this seems to close out a data block
        continue;
#count headers
    if line[0:4] == "****":
        header += 1;
        b_flag=1;
        continue;

# beginning NODE block
    if line[1:5] == "NODE":
        continue;

#beginning ELEMENT block
    if line[1:5] == "ELEM":
        continue;

#beginning of sideset block
    if line[1:6] == "ELSET":
        ss_count+=1;
        continue;

#beginning of nodeset block
    if line[1:5] == "NSET":
        ns_count += 1;
        ns_flag = 1;
        continue;

#Deal with end of sideset block
    if line[1:5] == "SURF":       
        continue;

    if line[0:3] == "SS1":       
        continue;
        
        
#Nodes
    if header == 2:
        lli = map(float,line.split(','));
        v_dict[lli[0]]=lli[1::];
        lli = lli[1::];
        if b_flag:
            print "processing vertices"
            vertex_array=numpy.array(lli,'f8');
            b_flag=0;
        else:
            vertex_array = numpy.vstack((vertex_array,lli));
        continue;

    if header == 3:
        lli = map(int,line.split(','));
        e_dict[lli[0]]=lli[1::];
#pflotran will need to know what kind of element 8 is hex
        lli[0] = 8;
        if b_flag:
            print "processing elements"
            element_array=numpy.array(lli,'i8');
            b_flag=0;
        else:
            element_array=numpy.vstack((element_array,lli));
        continue;
    
#nodset for a boundary condition
    if header == 4:
        if ns_count >= 1:
            if b_flag:
                nset = map(int,line.split(',')[0:-1]);
                b_flag=0;
            else:
                nset= nset+map(int,line.split(',')[0:-1]);
            continue;

#sideset (element set) for a boundary condition
    if header == 5:
        if ss_count >= 1:
            if b_flag:
                elset = map(int,line.split(',')[0:-1]);
                b_flag=0;
            else:
                elset= elset+map(int,line.split(',')[0:-1]);
            continue;

#loop through elements and write the face

if ns_flag:
	count = 1;
	nset = numpy.array(nset);
	for e in elset:
		nodes_e = numpy.array(e_dict[e]);
		face_e = numpy.array([],'int8');
		for n in nodes_e:
			face_e = numpy.concatenate((face_e,nset[numpy.where(nset==n)[0]]));
		
		if count == 1:
			faces = numpy.concatenate((numpy.array([4]),face_e));
		else:
			face_e = numpy.concatenate((numpy.array([4]),face_e));
			faces = numpy.vstack((faces,face_e));
		count += 1;
        
#pdb.set_trace();
h5dset = h5file.create_dataset('Domain/Cells',data=element_array);
h5dset = h5file.create_dataset('Domain/Vertices',data=vertex_array);
if ns_flag:
    h5dset = h5file.create_dataset('Regions/Tunnel_sideset',data=faces);

h5file.close();

