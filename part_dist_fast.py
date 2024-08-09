import os
import numpy as np
import pyvista as pv
import pickle
from tqdm import tqdm
from subdomain_functions import set_subdomain,set_part_subdomain,get_adjacent_indices,add_p_to_subdomains
import argparse

def rand_part(Lx,Ly,Lz,per_bc,rad):
    while True:
        part = np.random.rand(1,3)*np.array([Lx,Ly,Lz])
        if per_bc[0] == 0 and (part[0][0] < rad or part[0][0] > Lx-rad): continue
        if per_bc[1] == 0 and (part[0][1] < rad or part[0][1] > Ly-rad): continue
        if per_bc[2] == 0 and (part[0][2] < rad or part[0][2] > Lz-rad): continue
        break
    return part

def calc_dist_3periodic(part1, part2, Lx, Ly, Lz):
    offsets = np.array([
        [-Lx, -Ly, -Lz], [0, -Ly, -Lz], [Lx, -Ly, -Lz],
        [-Lx, 0, -Lz], [0, 0, -Lz], [Lx, 0, -Lz],
        [-Lx, Ly, -Lz], [0, Ly, -Lz], [Lx, Ly, -Lz],
        [-Lx, -Ly, 0], [0, -Ly, 0], [Lx, -Ly, 0],
        [-Lx, 0, 0], [0, 0, 0], [Lx, 0, 0],
        [-Lx, Ly, 0], [0, Ly, 0], [Lx, Ly, 0],
        [-Lx, -Ly, Lz], [0, -Ly, Lz], [Lx, -Ly, Lz],
        [-Lx, 0, Lz], [0, 0, Lz], [Lx, 0, Lz],
        [-Lx, Ly, Lz], [0, Ly, Lz], [Lx, Ly, Lz]
    ])

    part1_pos = part1 + offsets
    dists = np.linalg.norm(part1_pos - part2, axis=1)
    return np.min(dists)


#Initialize flags
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--radius', type=float, required=True,
                    help='Particles radii')
parser.add_argument('-d', '--delta', type=float, required=True,
                    help='Minimal distance between particles')
parser.add_argument('-p', '--phi', type=float, required=True,
                    help='Volume fraction of particles')
parser.add_argument('-l','--domain_length', type=str,required=True,
                    help='Size of the domain in the format of Lx:Ly:Lz')
parser.add_argument('-s','--subdomains_number', type=str,required=True,
                    help='Number of subdomains in the format of nx:ny:nz')
parser.add_argument('-per_bc', '--periodic_bc', type=str, required=True,
                    help='Periodic BC over axis, set in "0:1:1" format, where 1 correspond to periodic BC and\
                    0 to wall BC')
parser.add_argument('-o', '--out_folder', type=str, required=False,
                    help='Outpur folder for .pickle and .inp files')
parser.add_argument('-v', '--visualize', type=bool, required=False,
                    help='if true, provides a visualization by pyvista of resulting particles')
args = parser.parse_args()
    
    
#Input variables
Lx = float(args.domain_length.split(":")[0])
Ly = float(args.domain_length.split(":")[1])
Lz = float(args.domain_length.split(":")[2])

#Subdomain number
nx = int(args.subdomains_number.split(":")[0])
ny = int(args.subdomains_number.split(":")[1])
nz = int(args.subdomains_number.split(":")[2])

phi = args.phi
rad = args.radius
delta = args.delta

#Periodic BC
per_bc = [int(args.periodic_bc.split(":")[0]),
          int(args.periodic_bc.split(":")[1]),
          int(args.periodic_bc.split(":")[2])]

if args.out_folder:
    out_folder = args.out_folder + "/"
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)
        print("Output folder created")
    else: print("Output folder exists")
    
else: out_folder = ""
out_name_ptrn = f"p_mobile_phi_{str(phi)}_Lx_{str(Lx)}_Ly_{str(Ly)}_Lz_{str(Lz)}_r_{str(rad)}_delta_{str(delta)}"
out_inp_file =    out_folder + out_name_ptrn + ".inp"
out_pickle_file = out_folder + out_name_ptrn + ".pickle" 

#Calculate total number of particles
V_domain = Lx*Ly*Lz
n_part = int(phi*V_domain/(1.33333333*np.pi*rad**3))

part_list = []

#Add particles to the list
for i in tqdm(range(n_part),desc="Inserting particles"):
    while True:
        part = rand_part(Lx,Ly,Lz,per_bc,rad).flatten()
        break_flag = True

        if len(part_list) > 0:
            
            #Assign particles to corresponding subdomains
            if len(part_list) == 1:
                subdoms = set_subdomain(Lx,Ly,Lz,nx,ny,nz,part_list)
            
            #Find subdomain of master particle
            ii,jj,kk = set_part_subdomain(Lx, Ly, Lz, nx, ny, nz, part)
            adj_idx_list = get_adjacent_indices(ii, jj, kk, nx, ny, nz)
            
            
            for subdom_idx in adj_idx_list:
                iii,jjj,kkk = subdom_idx
                if len(subdoms[iii][jjj][kkk]) > 0:
                    for part_idx in subdoms[iii][jjj][kkk]:
                        pos = part_list[part_idx]
                        dist = calc_dist_3periodic(pos,part,Lx,Ly,Lz)
                        if dist < 2*(rad+delta): 
                            break_flag = False
                            break                     

        if break_flag == True: 
            part_list.append(part)
            if len(part_list) > 1: subdoms = add_p_to_subdomains(subdoms,Lx,Ly,Lz,nx,ny,nz,part,len(part_list)-1)
            break

with open(out_pickle_file, 'wb') as file:
    pickle.dump(part_list, file)

with open(out_inp_file,"w") as f:
    f.write(str(n_part)+"\n")
    for part in part_list:
        f.write(str(part[0])+"\t"+str(part[1])+"\t"+str(part[2])+"\t"+str(rad)+"\n")
        
if args.visualize:
    pl = pv.Plotter()
    domain = pv.Box(bounds=(0,Lx,0,Ly,0,Lz))
    pl.add_mesh(domain,style='wireframe')
    for part in part_list:
        sphere = pv.Sphere(center=part,radius=rad)
        pl.add_mesh(sphere)
    pl.show()
