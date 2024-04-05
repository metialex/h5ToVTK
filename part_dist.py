import numpy as np
import pyvista as pv
import pickle
from tqdm import tqdm
import argparse




parser = argparse.ArgumentParser()
parser.add_argument('-Lx', '--Lx', type=float, required=True, help='Lx')
parser.add_argument('-Ly', '--Ly', type=float, required=True, help='Ly')
parser.add_argument('-Lz', '--Lz', type=float, required=True, help='Lz')
parser.add_argument('-phi', '--phi', type=float, required=True, help='phi')
parser.add_argument('-rad', '--rad', type=float, required=True, help='rad')
parser.add_argument('-delta', '--delta', type=float, required=True, help='delta')
parser.add_argument('-per_bc', '--periodic_bc', type=str, required=True,
                    help='Periodic BC over axis, set in "1:1:1" format')
parser.add_argument('-xmin', '--xmin', type=float, required=False, help='xmin')
parser.add_argument('-ymin', '--ymin', type=float, required=False, help='ymin')
parser.add_argument('-zmin', '--zmin', type=float, required=False, help='zmin')
args = parser.parse_args()



#Input variables
Lx = args.Lx
Ly = args.Ly
Lz = args.Lz

phi = args.phi
rad = args.rad
delta = args.delta

per_bc = [int(args.periodic_bc.split(":")[0]),
          int(args.periodic_bc.split(":")[1]),
          int(args.periodic_bc.split(":")[2])]

out_file_name = "p_mobile.inp"
pickle_name = "p_mobile.pickle"

#Calculation
V_domain = Lx*Ly*Lz
n_part = int(phi*V_domain/(1.33333333*np.pi*rad**3))

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

part_list = []
for i in tqdm(range(n_part)):
    while True:
        part = rand_part(Lx,Ly,Lz,per_bc,rad)
        break_flag = True
        if len(part_list) > 0:
            for l_part in part_list:
                dist = calc_dist_3periodic(l_part,part,Lx,Ly,Lz)
                if dist < 2*(rad+delta): 
                    break_flag = False
                    break

        if break_flag == True: 
            part_list.append(part)
            break


pl = pv.Plotter()
domain = pv.Box(bounds=(0,Lx,0,Ly,0,Lz))
pl.add_mesh(domain,style='wireframe')
for part in part_list:
    sphere = pv.Sphere(center=part[0],radius=rad)
    pl.add_mesh(sphere)
pl.show()

#Translate positions of particles according to the domain
xmin = ymin = zmin = 0
if args.xmin: xmin = args.xmin
if args.xmin: xmin = args.xmin
if args.xmin: xmin = args.xmin

for part in part_list:
    part += [xmin,ymin,zmin]

with open(pickle_name, 'wb') as file:
    pickle.dump(part_list, file)

with open(out_file_name,"w") as f:
    f.write(str(n_part)+"\n")
    for part in part_list:
        f.write(str(part[0][0])+"\t"+str(part[0][1])+"\t"+str(part[0][2])+"\t"+str(rad)+"\n")
