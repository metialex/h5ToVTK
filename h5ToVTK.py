import numpy as np
import h5py
import glob
import os
import argparse
import struct
import sys
from tqdm import tqdm

def eulerH5toVTK(output_path, files):
        for input_file in files:
            tmp = input_file.split('Data_')[1]
            print('Processing - ' + tmp + ' file')
            file = output_path + '/Data_' + tmp.split('.h5')[0] + '.vtk'
            with h5py.File(input_file, 'r') as hdf:

                NX = np.array(hdf.get('grid').get('NX'))[0]
                NY = np.array(hdf.get('grid').get('NY'))[0]
                NZ = np.array(hdf.get('grid').get('NZ'))[0]

                dx = np.array(hdf.get('grid').get('xu'))
                dy = np.array(hdf.get('grid').get('yv'))
                dz = np.array(hdf.get('grid').get('zw'))

                pressure = np.array(hdf.get('p'))
                if(hdf.get('p_avg')):
                    pressure_av = np.array(hdf.get('p_avg'))
                else:
                    pressure_av = 0
                u = np.array(hdf.get('u'))
                v = np.array(hdf.get('v'))
                w = np.array(hdf.get('w'))

                try:
                    conc = np.array(hdf.get('Conc').get('0'))

                except:
                    conc = 0.0

 
                eulerVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w,pressure_av,conc)
def eulerVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w,pressure_av,conc):
    
    f = open(file,'wb')

    #General information about file
    f.write(b'# vtk DataFile Version 2.0\n')
    f.write(b'Sample rectilinear grid\n')
    f.write(b'BINARY\n')
    f.write(b'DATASET RECTILINEAR_GRID\n')

    #Grid characteristics
    tmp = 'DIMENSIONS ' + str(NX) + ' ' + str(NY) + ' ' + str(NZ) + '\n'
    sentence = bytearray(tmp.encode("ascii"))
    f.write(sentence)
    
    #X coordinates
    tmp = 'X_COORDINATES ' + str(NX) + ' float\n'
    sentence = bytearray(tmp.encode("ascii"))
    f.write(sentence)
    for i in dx:
        f.write(struct.pack(">f",i))
    f.write(b'\nMETADATA\nINFORMATION 0\n\n')
    
    #Y coordinates
    tmp = 'Y_COORDINATES ' + str(NY) + ' float\n'
    sentence = bytearray(tmp.encode("ascii"))
    f.write(sentence)
    for i in dy:
        f.write(struct.pack(">f",i))
    f.write(b'\nMETADATA\nINFORMATION 0\n\n')
    
    #Z coordinates
    tmp = 'Z_COORDINATES ' + str(NZ) + ' float\n'
    sentence = bytearray(tmp.encode("ascii"))
    f.write(sentence)
    for i in dz:
        f.write(struct.pack(">f",i))
    f.write(b'\nMETADATA\nINFORMATION 0\n\n')

    #Pressure field
    tmp = 'CELL_DATA ' + str((NX-1)*(NY-1)*(NZ-1)) + '\n'
    sentence = bytearray(tmp.encode("ascii"))
    f.write(sentence)
    f.write(b'SCALARS p_rk float\n')
    f.write(b'LOOKUP_TABLE default\n')
    for k in range(NZ-1):
        for j in range(NY-1):
            for i in range(NX-1):
                f.write(struct.pack(">f",pressure[k][j][i]))
    f.write(b'\nMETADATA\nINFORMATION 0\n\n')

    #concentration field
    try:
        print(conc[1][1][1])
        f.write(b'\n')
        f.write(b'SCALARS concentration float\n')
        f.write(b'LOOKUP_TABLE default\n')
        for k in range(NZ-1):
            for j in range(NY-1):
                for i in range(NX-1):
                    f.write(struct.pack(">f",conc[k][j][i]))
        f.write(b'\nMETADATA\nINFORMATION 0\n\n')
    except:
        print("No concentration")

    
    #Velocity field
    f.write(b'\n')
    f.write(b'VECTORS velocity float\n')
    for k in range(NZ-1):
        for j in range(NY-1):
            for i in range(NX-1):
                f.write(struct.pack(">f",u[k][j][i]))
                f.write(struct.pack(">f",v[k][j][i]))
                f.write(struct.pack(">f",w[k][j][i]))
    f.write(b'\nMETADATA\nINFORMATION 0\n')
    f.close()

def lagrangianH5toVTK(output_path, files):
        for input_file in files:
            tmp = input_file.split('Particle_')[1]
            print('Processing - ' + tmp + ' file')
            file = output_path + '/Particle_' + tmp.split('.h5')[0] + '.vtk'
            with h5py.File(input_file, 'r') as hdf:
                part_word = 'undefined'
                Position_1 = np.empty([1, 3])
                Position_2 = np.empty([1, 3])
                Radius_1 = np.empty([1, 1])
                Radius_2 = np.empty([1, 1])
                Velocity_1 = np.empty([1, 3])
                Velocity_2 = np.empty([1, 3])
                Omega_1 = np.empty([1, 3])
                Omega_2 = np.empty([1, 3])
                try:
                    
                    try:
                        list(hdf.get('mobile').items())
                        part_word = 'mobile'
                        Position_1 = np.array(hdf.get(part_word).get('X'))
                        Radius_1 = np.array(hdf.get(part_word).get('R'))
                        Velocity_1 = np.array(hdf.get(part_word).get('U'))
                        Omega_1 = np.array(hdf.get(part_word).get('Omega'))
                    except:
                        print('No mobile particles')
                    try:
                        list(hdf.get('fixed').items())
                        part_word = 'fixed'
                        Position_2 = np.array(hdf.get(part_word).get('X'))
                        Radius_2 = np.array(hdf.get(part_word).get('R'))
                        Velocity_2 = np.array(hdf.get(part_word).get('U'))
                        Omega_2 = np.array(hdf.get(part_word).get('Omega'))
                    except:
                        print('No fixed particles')
                    Position = np.concatenate((Position_1,Position_2))
                    Radius = np.concatenate((Radius_1,Radius_2))
                    Velocity = np.concatenate((Velocity_1,Velocity_2))
                    Omega = np.concatenate((Omega_1,Omega_2))
                    
                    lagrangianVTKwrite(file,Position, Radius, Velocity, Omega)
                except:
                    print('Lagrangian data is weird')         
def lagrangianVTKwrite(file,Position,Radius, Velocity, Omega):
    f = open(file,'w')

    #General information about file
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Sample rectilinear grid\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')


    #Position
    f.write('POINTS ' + str(len(Position)) + ' double\n')
    for i in Position:
        for j in i:
            f.write(str(j)+' ')
        f.write('\n')

    f.write('CELLS 0 0\nCELL_TYPES 0')
    #Radius
    f.write('POINT_DATA ' + str(len(Position)) + '\n')
    f.write('SCALARS R double\n')
    f.write('LOOKUP_TABLE default\n')
    for i in Radius:
        f.write(str(i[0]) + '\n')
    #Velocity
    f.write('\nVECTORS Velocity float\n')
    for i in Velocity:
        f.write(str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')
    #Angular velocity
    f.write('\nVECTORS Omega float\n')
    for i in Omega:
        f.write(str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')

    f.close()
def findFiles(path, keyword, index_min, index_max):
    line = path + keyword + '*' #keyword = '/Particle_* or '/Data''
    files = glob.glob(line)
    res = []
    for i in files:
        tmp = i.split(keyword)
        tmp = tmp[1].split('.h5')
        tmp = tmp[0]
        res.append(int(tmp))
        files_res =[]
    for i in range(len(files)):
        if res[i] <= index_max and res[i] >= index_min:
            files_res.append(files[i])
    files_res.sort()
    return files_res

def lagrangian_extract_quantity(output_path, quantity_list, files):
    #--------------------Determine dimensions of each quantity-------------------------#
    dimension_list = []
    tmp = files[0].split('Particle_')[1]
    with h5py.File(files[0], 'r') as hdf:
        part_word = 'undefined'
        try:
            list(hdf.get('mobile').items())
            part_word = 'mobile'
        except:
            try:
                list(hdf.get('fixed').items())
                part_word = 'fixed'
            except:
                print('Lagrangian Data is weird...')   
        for quantity in quantity_list:
            tmp1 = np.array(hdf.get(part_word).get(quantity))
            dimension_list.append(len(tmp1[0]))
    
        #--------------------Fill out the data array-------------------------#
    res = []#i - timestep, j- index of quantity, k -particle index
    time = []
    for input_file in files:
        tmp = input_file.split('Particle_')[1]
        print('Processing - ' + tmp + ' file')
        file = output_path + '/Particle_' + tmp.split('.h5')[0] + '.vtk'
        tmp2 = []
        with h5py.File(input_file, 'r') as hdf:
            time.append(hdf.get('time')[0])
            #print(time)
            part_word = 'undefined'
            try:
                list(hdf.get('mobile').items())
                part_word = 'mobile'
            except:
                try:
                    list(hdf.get('fixed').items())
                    part_word = 'fixed'
                except:
                    print('Lagrangian Data is weird...')
            
            #print(list(hdf.get('mobile').items()))
            tmp1 = []    
            for quantity in quantity_list:
                tmp = np.array(hdf.get(part_word).get(quantity))
                for i in range(len(tmp[0])):
                    for j in range(len(tmp)):
                        tmp1.append(tmp[j][i])
                    tmp2.append(tmp1)
                    tmp1 = [] 
        res.append(tmp2)

    #--------------------Write the array into the .dat file-------------------------#
    for q in range(len(quantity_list)):
        for qq in range(dimension_list[q]):
            file = quantity_list[q] + str(qq) + '.dat'
            f = open(file,'w')
            for i in range(len(res)):
                f.write(str(time[i]) + ',')
                for j in res[i][qq]:
                    f.write(str(j) + ',')
                f.write('\n')
    
def write_Reader_Vector_Velocity_xmf(path,no_files,Nx,Ny,Nz,time):
    with open("Reader_vector_python.xmf",'w') as f:
        Nx_i = Nx - 1
        Ny_i = Ny - 1
        Nz_i = Nz - 1
        N_vec = 3

        f.write(f'<?xml version=\"1.0\" ?>\n')
        f.write(f'<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n')
        f.write(f'<Xdmf Version=\"2.0\">\n')
        f.write(f'  <Domain>\n\n')

        f.write(f'    <Topology TopologyType="3DRectMesh" NumberOfElements="{Nz_i:d} {Ny_i:d} {Nx_i:d}"/>\n')
        f.write(f'    <Geometry GeometryType="VXVYVZ">\n')
        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nx_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Nx_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nx_i:d}">\n')
        f.write(f'          Vector_{time[0][0]:d}.h5:/grid/xc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Ny_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Ny_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Ny_i:d}">\n')
        f.write(f'          Vector_{time[0][0]:d}.h5:/grid/yc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Nz_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nz_i:d}">\n')
        f.write(f'          Vector_{time[0][0]:d}.h5:/grid/zc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n')
        f.write(f'    </Geometry>\n\n')

        f.write(f'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n\n')


        for i in range(no_files):
            f.write(f'      <Grid Name="SpatialGrid_{(i-1):d}" GridType="Uniform">\n', )
            f.write(f'        <Time Value="{time[i][1]:f}"/>\n')
            f.write(f'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
            f.write(f'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
            f.write(f'        <Attribute Name="vector" AttributeType="Vector" Center="Node">\n')
            f.write(f'          <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d} {Ny_i:d} {Nx_i:d} {N_vec:d}" Type="HyperSlab">\n')
            f.write(f'            <DataItem Dimensions="3 4" Format="XML">\n')
            f.write(f'              0    0    0    0 \n')
            f.write(f'              1    1    1    1 \n')
            f.write(f'              {Nz_i:d} {Ny_i:d} {Nx_i:d} {N_vec:d}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="{Nz_i:d} {Ny_i:d} {Nx_i:d} {N_vec:d}">\n')
            f.write(f'              Vector_{time[i][0]:f}.h5:/vector_velocity\n')
            f.write(f'            </DataItem>\n')
            f.write(f'          </DataItem>\n')
            f.write(f'        </Attribute>\n')
            f.write(f'      </Grid>\n\n')

        f.write(f'    </Grid>\n')
        f.write(f'  </Domain>\n')
        f.write(f'</Xdmf>\n')
        f.close()
        print('Wrote  Reader_vector.xmf\n')

def write_Reader_Scalar_Velocity_xmf(path, no_files, time, Nx, Ny, Nz, var):

    if var == "u":
        N = Nx
        Nx_i = Nx
        Ny_i = Ny - 1
        Nz_i = Nz - 1   
        x = 'xu'
        y = 'yc'
        z = 'zc'
    elif var == "v":
        N = Ny
        Nx_i = Nx - 1
        Ny_i = Ny
        Nz_i = Nz - 1 
        x = 'xc'
        y = 'yv'
        z = 'zc'
    elif var == "w":
        N = Nz
        Nx_i = Nx - 1
        Ny_i = Ny - 1
        Nz_i = Nz 
        x = 'xc'
        y = 'yc'
        z = 'zw'

    with open('Reader_'+var+'.xmf','w') as f:
    
        f.write(f'<?xml version=\"1.0\" ?>\n') 
        f.write(f'<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n')
        f.write(f'<Xdmf Version=\"2.0\">\n')
        f.write(f'  <Domain>\n\n')

        f.write(f'    <Topology TopologyType="3DRectMesh" NumberOfElements="{Nz_i:d} {Ny_i:d} {Nx_i:d}"/>\n')
        f.write(f'    <Geometry GeometryType="VXVYVZ">\n')
        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nx_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Nx_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nx:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/{x}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Ny_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Ny_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Ny:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/{y}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Nz_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nz_i:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/{z}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n')
        f.write(f'    </Geometry>\n\n')

        f.write(f'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n\n')


        for i in range(no_files):
            f.write(f'      <Grid Name="SpatialGrid_{i-1:d}" GridType="Uniform">\n')
            f.write(f'        <Time Value="{time[i][1]:f}"/>\n')
            f.write(f'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
            f.write(f'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
            f.write(f'        <Attribute Name="{var}" AttributeType="Scalar" Center="Node">\n')
            f.write(f'          <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d} {Ny_i:d} {Nx_i:d}" Type="HyperSlab">\n')
            f.write(f'            <DataItem Dimensions="3 3" Format="XML">\n')
            f.write(f'              0    0    0 \n')
            f.write(f'              1    1    1\n')
            f.write(f'              {Nz_i:d} {Ny_i:d} {Nx_i:d}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="{Nx:d} {Ny:d} {Nz:d}">\n')
            f.write(f'              Data_{time[i][0]:d}.h5:/{var}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'          </DataItem>\n')
            f.write(f'        </Attribute>\n')
            f.write(f'        <Attribute Name="vf{var}" AttributeType="Scalar" Center="Node">\n')
            f.write(f'          <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d} {Ny_i:d} {Nx_i:d}" Type="HyperSlab">\n')
            f.write(f'            <DataItem Dimensions="3 3" Format="XML">\n')
            f.write(f'              0    0    0 \n')
            f.write(f'              1    1    1\n')
            f.write(f'              {Nz_i:d} {Ny_i:d} {Nx_i:d}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="{Nx:d} {Ny:d} {Nz:d}">\n')
            f.write(f'              Data_{time[i][0]:d}.h5:/vf{var}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'          </DataItem>\n')
            f.write(f'        </Attribute>\n')
            f.write(f'      </Grid>\n\n')


        f.write(f'    </Grid>\n')
        f.write(f'  </Domain>\n')
        f.write(f'</Xdmf>\n')

    print('Wrote  Reader_'+var+'.xmf\n')

def write_Reader_Particle_xmf(path, no_files, time, type, Type, Np, att_list):

    with open('Reader_p_'+type+'.xmf','w') as f:

        f.write(f'<?xml version="1.0" ?>\n')
        f.write(f'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write(f'<Xdmf Version="2.0">\n')
        f.write(f'  <Domain>\n\n')
        f.write(f'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n')

        for i in range(no_files):
            f.write(f'      <Grid Name="{Type}Grid_{i-1:d}">\n')
            f.write(f'        <Time Value="{time[i][1]:f}"/>\n')
            f.write(f'        <Topology Type="Polyvertex" NumberOfElements="{Np:d}" />\n\n') 
            f.write(f'        <Geometry Type="XYZ">\n')
            f.write(f'          <DataItem Format="HDF" Dimensions="{Np:d} 3">\n')
            f.write(f'            Particle_{time[i][0]:d}.h5:/{type}/X\n')
            f.write(f'          </DataItem>\n')
            f.write(f'        </Geometry>\n')
            
            for j in range(len(att_list)):
                f.write(f'\n        <Attribute Name="{att_list[j][0]}" AttributeType="{att_list[j][2]}" Center="Node">\n')
                f.write(f'          <DataItem Format="HDF" NumberType="Double" Dimensions="{Np:d} {att_list[j][1]:d}">\n')
                f.write(f'            Particle_{time[i][0]:d}.h5:/{type}/{att_list[j][0]}\n')
                f.write(f'          </DataItem>\n')
                f.write(f'        </Attribute>\n')

            f.write(f'      </Grid>\n\n')


        f.write(f'    </Grid>\n')
        f.write(f'  </Domain>\n')
        f.write(f'</Xdmf>\n')

    print('Wrote  Reader_p_'+type+'.xmf\n')

def write_Reader_scalar_xmf(path, no_files, time, Nx, Ny, Nz, var,var_name):

    Nx_i = Nx - 1
    Ny_i = Ny - 1
    Nz_i = Nz - 1   


    with open(f'Reader_{var_name}.xmf','w') as f:

        f.write(f'<?xml version=\"1.0\" ?>\n') 
        f.write(f'<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n')
        f.write(f'<Xdmf Version=\"2.0\">\n')
        f.write(f'  <Domain>\n\n')

        f.write(f'    <Topology TopologyType="3DRectMesh" NumberOfElements="{Nz_i:d} {Ny_i:d} {Nx_i:d}"/>\n')
        f.write(f'    <Geometry GeometryType="VXVYVZ">\n')
        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nx_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1    {Nx_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nx:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/xc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Ny_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1   {Ny_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Ny:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/yc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n\n')

        f.write(f'      <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d}" Type="HyperSlab">\n')
        f.write(f'        <DataItem Dimensions="3" Format="XML">\n')
        f.write(f'          0    1   {Nz_i:d}\n')
        f.write(f'        </DataItem>\n')
        f.write(f'        <DataItem Format="HDF" Dimensions="{Nz:d}">\n')
        f.write(f'          Data_{time[0][0]:d}.h5:/grid/zc\n')
        f.write(f'        </DataItem>\n')
        f.write(f'      </DataItem>\n')
        f.write(f'    </Geometry>\n\n')

        f.write(f'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n\n')

        for i in range(no_files):

            f.write(f'      <Grid Name="SpatialGrid_{i-1:d}" GridType="Uniform">\n')
            f.write(f'        <Time Value="{time[i][1]:f}"/>\n')
            f.write(f'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
            f.write(f'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
            f.write(f'        <Attribute Name="{var_name}" AttributeType="Scalar" Center="Node">\n')
            f.write(f'          <DataItem ItemType="HyperSlab" Dimensions="{Nz_i:d} {Ny_i:d} {Nx_i:d}" Type="HyperSlab">\n')
            f.write(f'            <DataItem Dimensions="3 3" Format="XML">\n')
            f.write(f'              0    0    0 \n')
            f.write(f'              1    1    1\n')
            f.write(f'              {Nz_i:d} {Ny_i:d} {Nx_i:d}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="{Nx:d} {Ny:d} {Nz:d}">\n')
            f.write(f'              Data_{time[i][0]:d}.h5:/{var}\n')
            f.write(f'            </DataItem>\n')
            f.write(f'          </DataItem>\n')
            f.write(f'        </Attribute>\n')
            f.write(f'      </Grid>\n\n')


        f.write(f'    </Grid>\n')
        f.write(f'  </Domain>\n')
        f.write(f'</Xdmf>\n')



    print(f'Wrote  Reader_{var_name}.xmf\n')

def if_exists(file,subdir):
    with h5py.File(file, 'r') as hdf:
        if hdf.get(subdir) is not None:
            return 1
        else:
            return 0

def attributes(file,subdir):
    result = []
    with h5py.File(file, 'r') as hdf:
        attributes = list(hdf.get(subdir))
        for attribute in attributes:
            size = len(hdf.get(subdir).get(attribute)[0])
            if size == 3: type = 'Vector'
            elif size == 1: type = 'Scalar'
            elif size == 9: type = 'Tensor'
            else: continue
            result.append([attribute,size,type])
    return result

def numper_of_particles(file,subdir):
    with h5py.File(file, 'r') as hdf:
        attributes = list(hdf.get(subdir))
        Np = len(hdf.get(subdir).get(attributes[0]))
    return Np

def set_time(files):
    res = []
    for file in files:
        tmp = file.split('_')[1]
        idx = int(tmp.split('.h5')[0])
        with h5py.File(file, 'r') as hdf:
            time = hdf.get('time')[0]
        res.append([idx,time])
    return res

def set_Nx_Ny_Nz(file):
    with h5py.File(file, 'r') as hdf:
        NX = hdf.get('grid').get('NX')[0]
        NY = hdf.get('grid').get('NY')[0]
        NZ = hdf.get('grid').get('NZ')[0]
    return NX,NY,NZ

def write_Rotn_paraview(files_l):
    for file_name in tqdm(files_l):
        with h5py.File(file_name,"a") as f:
            #delete variable if exists
            if "mobile/Rotn_paraview" in f:
                del f["mobile/Rotn_paraview"]
            
            Rotn_paraview = f["mobile/Rotn"][:]
            
            for jj in range(len(f["mobile/Rotn"])):
                Rotn_paraview[jj][1] = f["mobile/Rotn"][jj][3]
                Rotn_paraview[jj][2] = f["mobile/Rotn"][jj][6]
                Rotn_paraview[jj][3] = f["mobile/Rotn"][jj][1]
                Rotn_paraview[jj][5] = f["mobile/Rotn"][jj][7]
                Rotn_paraview[jj][6] = f["mobile/Rotn"][jj][2]
                Rotn_paraview[jj][7] = f["mobile/Rotn"][jj][5]
                
            f.create_dataset("mobile/Rotn_paraview", data=Rotn_paraview)


parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group()
group.add_argument('-v','--vtk',action='store_true',help='Convert .h5 data into .vtk files')
group.add_argument('-x','--xdmf',action='store_true', help='Create .xmf files for .h5 files')
parser.add_argument('-i', '--index', type=str,required=False, help='index range of data to convert (e.g. -i 2:5)')
parser.add_argument('-d', '--dat', type=str,required=False, help='Converting .h5 file into .dat file.\
     Avaliable variables - F_IBM, F_coll, F_rigid, Fc, Int_Omega_old,\
          Int_U_old, Omega, R, Rotn, T_IBM, T_coll, T_rigid, Tc, U, X, collision_data')
#group = parser.add_mutually_exclusive_group()
parser.add_argument('-el', '--eulerLagrangian', action='store_true', help='Converting Eulerian and Lagrangian data')
parser.add_argument('-e', '--euler', action='store_true', help='Converting only Eulerian data')
parser.add_argument('-l', '--lagrangian', action='store_true', help='Converting only Lagrangian data')
parser.add_argument('-rotn','--rotn_paraview',action='store_true',
                   help='Create additional Rotn tensor for visualizing rotation in paraview')
parser.add_argument('-as', '--add_scalar', type=str,required=False, help='')

args = parser.parse_args()

#defaul paths
input_path = ''
output_path = 'vtk_out/'

#defaul indexes
index_min = 0
index_max = 1e+10

#Keywords
euler_keyword = 'Data_'
lagrang_keyword = 'Particle_'

if args.xdmf:
    output_path = ''
    if args.index:
            index_min = int(args.index.split(':')[0])
            index_max = int(args.index.split(':')[1])
    #Create Lagrangian dataset
    try:
        files_l = findFiles(input_path,lagrang_keyword, index_min, index_max)

        #The information about number of particles and attributes is taken from the 1st Particle.h5 file
        #mobile particles
        if args.rotn_paraview:
            print('Writing R_otn_paraview')
            write_Rotn_paraview(files_l)
        
        if if_exists(files_l[0],'mobile'):
            att_list = attributes(files_l[0],'mobile')
            Np = numper_of_particles(files_l[0],'mobile')
            time = set_time(files_l)
            write_Reader_Particle_xmf(output_path, len(files_l), time, 'mobile', 'Mobile', Np, att_list)


        #fixed particles
        if if_exists(files_l[0],'fixed'):
            att_list = attributes(files_l[0],'fixed')
            Np = numper_of_particles(files_l[0],'fixed')
            time = set_time(files_l)
            write_Reader_Particle_xmf(output_path, len(files_l), time, 'fixed', 'Fixed', Np, att_list)
 
    except:
        print("Lagrangian files are not found")
    
    try:
        #Eulerian Fields
        files_e = findFiles(input_path,euler_keyword, index_min, index_max)
    except:
        print("Eulerian files are not found")
        sys.exit()
        
    time = set_time(files_e)
    Nx,Ny,Nz = set_Nx_Ny_Nz(files_e[0])
    #Velocity field
    write_Reader_Scalar_Velocity_xmf(output_path, len(files_e), time, Nx, Ny, Nz, 'u')
    write_Reader_Scalar_Velocity_xmf(output_path, len(files_e), time, Nx, Ny, Nz, 'v')
    write_Reader_Scalar_Velocity_xmf(output_path, len(files_e), time, Nx, Ny, Nz, 'w')
    #Pressure field
    write_Reader_scalar_xmf(output_path, len(files_e), time, Nx, Ny, Nz, 'p','pressure')
    if args.add_scalar is not None:
        if if_exists(files_e[0],args.add_scalar):
            write_Reader_scalar_xmf(output_path, len(files_e), time, Nx, Ny, Nz, args.add_scalar, 'C')
    
elif args.vtk:
    #Create output directory
    directory = 'vtk_out/'
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    if args.dat:
        quantity_list = []
        for i in (args.dat.split(',')):
            quantity_list.append(i)
        files = findFiles(input_path,lagrang_keyword, index_min, index_max)
        lagrangian_extract_quantity(directory, quantity_list, files)
        
    elif args.eulerLagrangian:
        if args.index:
            index_min = int(args.index.split(':')[0])
            index_max = int(args.index.split(':')[1])

        try:
            print('Eulerian data start')
            files_e = findFiles(input_path,euler_keyword, index_min, index_max)
            eulerH5toVTK(output_path, files_e)
            print('Eulerian data stop')
        except:
            print ('No Eulerian data')
        try:
            print('Lagrangian data start')
            files_l = findFiles(input_path,lagrang_keyword, index_min, index_max)
            lagrangianH5toVTK(output_path, files_l)
            print('Lagrangian data stop')
        except:
            print ('No Lagrangian data')
    elif args.euler:
        if args.index:
            index_min = int(args.index.split(':')[0])
            index_max = int(args.index.split(':')[1])

        try:
            print('Eulerian data start')
            files_e = findFiles(input_path,euler_keyword, index_min, index_max)
            eulerH5toVTK(output_path, files_e)
            print('Eulerian data stop')
        except:
            print ('No Eulerian data')
    elif args.lagrangian:
        if args.index:
            index_min = int(args.index.split(':')[0])
            index_max = int(args.index.split(':')[1])
        try:
            print('Lagrangian data start')
            files_l = findFiles(input_path,lagrang_keyword, index_min, index_max)
            lagrangianH5toVTK(output_path, files_l)
            print('Lagrangian data stop')
        except:
            print ('No Lagrangian data')
    else:
        if args.index:
            index_min = int(args.index.split(':')[0])
            index_max = int(args.index.split(':')[1])

        try:
            print('Eulerian data start')
            files_e = findFiles(input_path,euler_keyword, index_min, index_max)
            eulerH5toVTK(output_path, files_e)
            print('Eulerian data stop')
        except:
            print ('No Eulerian data')
        try:
            print('Lagrangian data start')
            files_l = findFiles(input_path,lagrang_keyword, index_min, index_max)
            lagrangianH5toVTK(output_path, files_l)
            print('Lagrangian data stop')
        except:
            print ('No Lagrangian data')

