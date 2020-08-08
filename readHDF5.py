import numpy as np
import h5py
#import tkinter
import glob


def rectilinearVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w):
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
    for i in range(NX):
        f.write('\n')
        for j in range(NY):
            f.write('\n')
            for k in range(NZ):
                f.write(str(pressure[i][j][k]) + ' ')
    #Velocity field
    f.write('\n')
    f.write('VECTORS velocity float\n')
    for i in range(NX):
        f.write('\n')
        for j in range(NY):
            f.write('\n')
            for k in range(NZ):
                f.write(str(u[i][j][k]) + ' ')
                f.write(str(v[i][j][k]) + ' ')
                f.write(str(w[i][j][k]) + ' ')
    f.close()
def eulerH5toVTK(path, files):
        for input_file in files:
            tmp = input_file.split('Data_')[1]
            print('Processing - ' + tmp + ' file')
            file = 'output/Data_' + tmp.split('.h5')[0] + '.vtk'
            with h5py.File(input_file) as hdf:
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
 
                rectilinearVTKwrite(file,NX,NY,NZ,dx,dy,dz,pressure,u,v,w)
def findFiles(path):
    line = path + '/Data_*'
    files = glob.glob(line)
    print(files)
    return files
    

path = '/home/metialex/Desktop/Work/PARTIES/Simulation_cases/myTestSimulation'
files = findFiles(path)
eulerH5toVTK(path, files)

#GUI 
#window = tkinter.Tk()
#window.geometry("1000x600")
#window.title('.h5 to .vtk')

#label
#l1 = tkinter.Label(window, text='Directory path')
#l1.grid(column=0,row=0)
#l1.place()

#Entry for folder path
#txt=tkinter.Entry(window,width=40)
#txt.grid(column=1,row=0)

#button
#bt1 = tkinter.Button(window, text="run")
#bt1.grid(column=0,row=1)
#window.mainloop()

