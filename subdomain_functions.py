import numpy as np

def set_part_subdomain(Lx:float,
                       Ly:float,
                       Lz:float,
                       nx:int,
                       ny:int,
                       nz:int,
                       p):
    """
    Sets the 3d index of subdomain array
    where particle p is located.
    
    Args:
        Lx,Ly,Lz (float)  - size of domain
        nx,ny,nz (int) - number of subdomain divisions
        p (list or np.array) - particle position vector 
    
    Returns:
        i,j,k (int) - indexes for appropriate subdomain where particle p is located
    """
        
    i = int(p[0]/Lx*nx)
    j = int(p[1]/Ly*ny)
    k = int(p[2]/Lz*nz)
    
    return i,j,k


def set_subdomain(Lx:float,
                  Ly:float,
                  Lz:float,
                  nx:int,
                  ny:int,
                  nz:int,
                  p_list):
    """
    Sets the subdomain 3D array, where particles indexes are
    distributed in corresponding elements of the subdomain
    
    Args:
        Lx,Ly,Lz (float): size of domain
        nx, ny, nz (int): number of subdomain divisions
        p_list (_type_): list of particle positions

    Returns:
        _type_: _description_
    """
    # Initialize the 3D list with zeros
    subdoms = [[[[] for k in range(nz)] for j in range(ny)] for i in range(nx)]
    
    #Distribute particles to appropriate subdomains
    for idx,p in enumerate(p_list):
        i,j,k = set_part_subdomain(Lx,Ly,Lz,nx,ny,nz,p)
        subdoms[i][j][k].append(idx)
    
    return subdoms

def add_p_to_subdomains(subdoms,
                        Lx:float,
                        Ly:float,
                        Lz:float,
                        nx:int,
                        ny:int,
                        nz:int,
                        p,
                        p_idx:int):
    """_summary_

    Args:
        subdoms (_type_): _description_
        Lx,Ly,Lz (float): size of domain
        nx, ny, nz (int): number of subdomain divisions
        p (list of np.array): particle position vector
        p_idx (int): particle index

    Returns:
        list: 3D list of subdomains where which element contain
        a list of particle indexes corresponding to subdomain
    """
    i,j,k = set_part_subdomain(Lx,Ly,Lz,nx,ny,nz,p)
    subdoms[i][j][k].append(p_idx)
    return subdoms

def get_adjacent_indices(i, j, k, nx, ny, nz):
    """
    Determine the neighboring indices in a 3D grid while considering triple periodic boundary conditions.

    Parameters:
    i, j, k (int): The current indices in the 3D grid for which adjacent indices are to be found.
    nx, ny, nz (int): The maximum indices in the x, y, and z dimensions of the 3D grid, respectively.

    Returns:
    list of lists: A list of lists where each sublist contains three integers representing the indices of an adjacent point in the 3D grid, considering periodic boundaries.

    Functionality:
    1. Periodic Boundary Calculation:
       - Uses a helper function, periodic_boundary, which calculates the previous, current, and next indices for a given index,
         ensuring that they wrap around the grid using modulo operations. This creates the effect of a toroidal structure where the grid edges are connected.
    2. Generating Neighboring Indices:
       - Computes the neighbors for each dimension (i, j, k) using the periodic_boundary function.
       - Generates all combinations of these neighboring indices using nested loops and returns them as a list.
    """
    
    # Helper function to calculate the periodic boundary
    def periodic_boundary(idx, max_idx):
        return (idx - 1) % max_idx, idx % max_idx, (idx + 1) % max_idx
    
    # Get the periodic boundary indices for each dimension
    i_neighbors = periodic_boundary(i, nx)
    j_neighbors = periodic_boundary(j, ny)
    k_neighbors = periodic_boundary(k, nz)
    
    # Generate all combinations of the neighboring indices
    res = []
    for ii in i_neighbors:
        for jj in j_neighbors:
            for kk in k_neighbors:
                res.append([ii, jj, kk])
    
    return res

    
    