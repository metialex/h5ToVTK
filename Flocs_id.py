import numpy as np
import h5py
import pandas as pd
import csv
import pyvista as pv
import glob
import math
import os
from tqdm import tqdm
import argparse


def compute_adjacency_list(df, cutoff_distance, Lx=None, Ly=None, Lz=None):
    
    # Initialize adjacency list
    adjacency_list = {}
    
    # Compute adjacency list
    pos = df[["x", "y", "z"]].to_numpy()
    rad = df["r"].to_numpy()
    for i, part in enumerate(pos):
        dist = pos - part
        #Check how exactly periodic BC are considered
        if Lx:
            dist[:, 0] -= np.round(dist[:, 0] / Lx) * Lx
        if Ly:
            dist[:, 1] -= np.round(dist[:, 1] / Ly) * Ly
        if Lz:
            dist[:, 2] -= np.round(dist[:, 2] / Lz) * Lz
        dist = (dist * dist).sum(axis=1)
        cutoff = rad[i] + rad + cutoff_distance
        is_adj = dist < cutoff * cutoff
        is_adj[i] = False
        adjacency_list[i] = np.where(is_adj)[0]
    return adjacency_list
    
def find_connected_groups_modified(adjacency_list):
    visited = set()
    connected_groups = []
    for start_node in adjacency_list:
        if start_node not in visited:
            connected_group = []
            stack = [start_node]
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    connected_group.append(node)
                    stack.extend(neighbor for neighbor in adjacency_list[node] if neighbor not in visited)
            connected_groups.append(connected_group)
    return connected_groups

def write_floc_data(part_h5_file_name, Lx,Ly,Lz,cutoff_distance):
    #Open h5 particle file and 
    with h5py.File(part_h5_file_name,'r+') as f:
        part_data = pd.DataFrame(f["mobile/X"][:], columns=["x", "y", "z"])
        part_data["r"] = f["mobile/R"][:]
    part_data["floc_size"] = 0
    part_data["floc_id"] = 0

    adjacency_list = compute_adjacency_list(part_data, cutoff_distance, Lx=Lx, Ly=Ly, Lz=Lz)
    connected_groups = find_connected_groups_modified(adjacency_list)

    #Set the size of the floc to every particle
    for idx,group in enumerate(connected_groups):
        for part in group:
            part_data.loc[part,"floc_size"] = [len(group)]
            part_data.loc[part,"floc_id"] = [idx]

    with h5py.File(part_h5_file_name,'a') as f:
        #Writing floc size
        if "mobile/floc_size" in f:
            floc_arr = part_data["floc_size"].values.reshape((len(part_data),1))
            f["mobile/floc_size"][:] = floc_arr
        else:
            floc_arr = part_data["floc_size"].values.reshape((len(part_data),1))
            f.create_dataset("mobile/floc_size", data=floc_arr)
        #Writing floc id
        if "mobile/floc_id" in f:
            floc_arr = part_data["floc_id"].values.reshape((len(part_data),1))
            f["mobile/floc_id"][:] = floc_arr
        else:
            floc_arr = part_data["floc_id"].values.reshape((len(part_data),1))
            f.create_dataset("mobile/floc_id", data=floc_arr)

def extract_domain_from_data(particle_data):
    with h5py.File(particle_data,"r+") as f:
        xmin = f["domain/xmin"][0]
        xmax = f["domain/xmax"][0]
        ymin = f["domain/ymin"][0]
        ymax = f["domain/ymax"][0]
        zmin = f["domain/zmin"][0]
        zmax = f["domain/zmax"][0]
        bc = f["domain/periodic"][:]
        return xmin,xmax,ymin,ymax,zmin,zmax,bc
    
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--cutoff_distance',
                    type=float, required=False,
                    help='Particles will be considered as part of the same floc, if the distance between their centers is smaller than Ri+Rj+coutoff_distance')
parser.add_argument('-i','--particle_index_range', type=str,required=False,help='Certain indexes of Particles_.h5 files to be analysed. The format if input is following -i 0:10')
args = parser.parse_args()

particle_data = glob.glob("Particle*")    
xmin,xmax,ymin,ymax,zmin,zmax,bc = extract_domain_from_data(particle_data[0])

Lx = Ly = Lz = False

if bc[0] == 1: Lx = xmax - xmin
if bc[1] == 1: Ly = ymax - ymin
if bc[2] == 1: Lz = zmax - zmin

#Default value of cutoff distance
cutoff_dist = 0.05
if args.cutoff_distance:
    cutoff_dist = args.cutoff_distance
if args.particle_index_range:
    i_min = int(args.particle_index_range.split(":")[0])
    i_max = int(args.particle_index_range.split(":")[1])
    particle_data = []
    for i in range(i_min,i_max,1):
        particle_data.append(f"Particle_{i}.h5")
else:
    particle_data = glob.glob("Particle*")
print(cutoff_dist)

for file in tqdm(particle_data):
    write_floc_data(file,Lx,Ly,Lz,cutoff_dist)
