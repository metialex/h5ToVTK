import numpy as np
import h5py
import sys
import glob
import os

def eulerH5toVTK(output_path, files):
        for input_file in files:
            tmp = input_file.split('Data_')[1]
            print('Processing - ' + tmp + ' file')
            file = output_path + '/Data_' + tmp.split('.h5')[0] + '.vtk'
            with h5py.File(input_file, 'r') as hdf:
                #ls = list(hdf.items())
                #print('List of datasets: \n',ls)
                NX = np.array(hdf.get('grid').get('NX'))[0]
                NY = np.array(hdf.get('grid').get('NY'))[0]
                NZ = np.array(hdf.get('grid').get('NZ'))[0]

                dx = np.array(hdf.get('grid').get('xc'))
                dy = np.array(hdf.get('grid').get('yc'))
                dz = np.array(hdf.get('grid').get('zc'))

                pressure = np.array(hdf.get('p'))
                u = np.array(hdf.get('u'))
                v = np.array(hdf.get('v'))
                w = np.array(hdf.get('w'))
 
                eulerVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w)
def eulerVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w):
    f = open(file,'w')

    #General information about file
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Sample rectilinear grid\n')
    f.write('ASCII\n')
    f.write('DATASET RECTILINEAR_GRID\n')

    #Grid characteristics
    f.write('DIMENSIONS ' + str(NX+1) + ' ' + str(NY+1) + ' ' + str(NZ+1) + '\n')
    #X coordinates
    f.write('X_COORDINATES ' + str(NX+1) + ' float\n')
    f.write('0.0 ')
    for i in dx:
        f.write(str(i)+' ')
    f.write('\n')
    #Y coordinates
    f.write('Y_COORDINATES ' + str(NY+1) + ' float\n')
    f.write('0.0')
    for i in dy:
        f.write(str(i)+' ')
    f.write('\n')
    #Z coordinates
    f.write('Z_COORDINATES ' + str(NZ+1) + ' float\n')
    f.write('0.0 ')
    for i in dz:
        f.write(str(i)+' ')
    f.write('\n')

    #Pressure field
    f.write('CELL_DATA ' + str(NX*NY*NZ) + '\n')
    f.write('SCALARS pressure float\n')
    f.write('LOOKUP_TABLE default\n')
    for k in range(NZ):
        f.write('\n')
        for j in range(NY):
            f.write('\n')
            for i in range(NX):
                f.write(str(pressure[k][j][i]) + ' ')
    #Velocity field
    f.write('\n')
    f.write('VECTORS velocity float\n')
    for k in range(NZ):
        f.write('\n')
        for j in range(NY):
            f.write('\n')
            for i in range(NX):
                f.write(str(u[k][j][i]) + ' ')
                f.write(str(v[k][j][i]) + ' ')
                f.write(str(w[k][j][i]) + ' ')
    f.close()

def lagrangianH5toVTK(output_path, files):
        for input_file in files:
            tmp = input_file.split('Particle_')[1]
            print('Processing - ' + tmp + ' file')
            file = output_path + '/Particle_' + tmp.split('.h5')[0] + '.vtk'
            with h5py.File(input_file, 'r') as hdf:
                ls = list(hdf.get('mobile').items())
                #print('List of datasets: \n',ls)
                Position = np.array(hdf.get('mobile').get('X'))
                Radius = np.array(hdf.get('mobile').get('R'))
                
                Velocity = np.array(hdf.get('mobile').get('U'))
                Omega = np.array(hdf.get('mobile').get('Omega'))
                #print(Velocity)
                lagrangianVTKwrite(file,Position, Radius, Velocity, Omega)          
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
    line = path + keyword + '*' #keyword = '/Particle_*'
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
    return files_res

    
with h5py.File('Particle_1332.h5', 'r') as hdf:
    ls = list(hdf.items())
    print('List of datasets: \n',ls)
    ls = list(hdf.get('mobile').items())
    for i in ls:
        print(i)

'''
print(sys.argv)

#Create output directory
directory = 'vtk_out/'
try:
    os.stat(directory)
except:
    os.mkdir(directory)

#initializing default variables
euler_keyword = 'Data_'
lagrang_keyword = 'Particle_'

input_path = ''#'/home/vowinckel/Desktop/PhD/Simulation_cases/72particles/Case1'
output_path = 'vtk_out/'#'/home/vowinckel/Desktop/PhD/Simulation_cases/72particles/Case1/vtk_out'

index_min = 0
index_max = 1e+10

#reading initial arguments
if len(sys.argv) == 1:
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
elif sys.argv[2] == '-index':
    index_min = int(sys.argv[3].split(':')[0])
    index_max = int(sys.argv[3].split(':')[1])

    if sys.argv[1] == '-e' or sys.argv[1] == '-el':  
        try:
            print('Eulerian data start')
            files_e = findFiles(input_path,euler_keyword, index_min, index_max)
            eulerH5toVTK(output_path, files_e)
            print('Eulerian data stop')
        except:
            print ('No Eulerian data')
    if sys.argv[1] == '-l' or sys.argv[1] == '-el':         
        try:
            print('Lagrangian data start')
            files_l = findFiles(input_path,lagrang_keyword, index_min, index_max)
            lagrangianH5toVTK(output_path, files_l)
            
        except:
            print ('No Lagrangian data')
else:
    print('Error in initial arguments\n')
    print('The example of initial arguments:')
    print('-el -index startIndex:stopIndex')
'''



