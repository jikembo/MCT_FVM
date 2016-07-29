import numpy as np
import Jacobian as jp

def coordinate(fa):
    """
	Read in Coordinate
    return npt: number pf points
           coor: coordinates 
	"""
    f_coor = open(fa)
#    npt = np.genfromtxt(f_coor, dtype = 'int', usecols = 0, max_rows = 1)
    coor = np.genfromtxt(f_coor, dtype = 'float', usecols =(1,2,3), skip_header = 0)
    npt = coor.shape[0]
    f_coor.close()
    return npt, coor
    
##########################################################################################################
    
def boundary_connectivity(fb, fc):
    """
    Read in Boundary and Connectivity
    return nbound: number of the boundary conditions
           nele: number of volume elements
           bound: connectivity for boundaries
           ijk: connectivity of volume elements
    """
    f_bound = open(fb)
    f_connect = open(fc)
#   ncond = np.genfromtxt(f_bound, dtype = 'int', usecols = 0, max_rows = 1)
    bound = np.genfromtxt(f_bound, dtype = 'int', usecols =(4,5,6,7,8), skip_header = 0)
    nbound = bound.shape[0]
#   nele = ncond - nbound
    ijk = np.genfromtxt(f_connect, dtype = 'int', usecols =(5,6,7,8,9,10,11,12))
    nele = ijk.shape[0]
    f_bound.close()
    f_connect.close()
    return nbound, nele, bound, ijk
    
