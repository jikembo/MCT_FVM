import numpy as np
def print_restart(Q_res, n_el, Re):
    f_restart = open(Re, 'w')
    for i in range(0,n_el,1):
        print >> f_restart, Q_res[i,0], Q_res[i,1], Q_res[i,2], Q_res[i,3], Q_res[i,4], Q_res[i,5], Q_res[i,6], Q_res[i,7]
    f_restart.close()

def print_output(Q, n_ele, inertia, k, i_step, coord, n_pt, ijk):
    str1 = "./Results/T"
    str2 = str(i_step)
    str3 = ".vtk"
    step_header = "".join((str1,str2,str3))
    f_out = open (step_header, 'w') 
    print >> f_out, "# vtk DataFile Version 2.0"
    print >> f_out, str2, "-th Step"
    print >> f_out, "ASCII"
    print >> f_out, "DATASET UNSTRUCTURED_GRID"
    print >> f_out, "POINTS", n_pt, "float"
    for i in range (0,n_pt,1):
        print >> f_out, coord[i,0], coord[i,1], coord[i,2]
    print >> f_out, "CELLS", n_ele, 9 * n_ele
    for i in range (0,n_ele,1):
        print >> f_out, "8", ijk[i,0]-1, ijk[i,1]-1, ijk[i,5]-1, ijk[i,4]-1, ijk[i,3]-1, ijk[i,2]-1, ijk[i,6]-1, ijk[i,7]-1
    print >> f_out, "CELL_TYPES", n_ele
    for i in range (0,n_ele,1):
        print >> f_out, "12"
    rho = np.zeros(n_ele, dtype=float)
    u = np.zeros(n_ele, dtype=float)
    v = np.zeros(n_ele, dtype=float)
    w = np.zeros(n_ele, dtype=float)
    ox = np.zeros(n_ele, dtype=float)
    oy = np.zeros(n_ele, dtype=float)
    oz = np.zeros(n_ele, dtype=float)
    temp = np.zeros(n_ele, dtype=float)
    ek = np.zeros(n_ele, dtype=float)
    p = np.zeros(n_ele, dtype=float)
    rho[0:n_ele] = Q[0:n_ele, 0]
    u[0:n_ele] = Q[0:n_ele, 1] / Q[0:n_ele, 0]
    v[0:n_ele] = Q[0:n_ele, 2] / Q[0:n_ele, 0]
    w[0:n_ele] = Q[0:n_ele, 3] / Q[0:n_ele, 0]
    ox[0:n_ele] = Q[0:n_ele, 4] / Q[0:n_ele, 0] / inertia
    oy[0:n_ele] = Q[0:n_ele, 5] / Q[0:n_ele, 0] / inertia
    oz[0:n_ele] = Q[0:n_ele, 6] / Q[0:n_ele, 0] / inertia
    ek[0:n_ele] = ((u[0:n_ele] * u[0:n_ele]) + (v[0:n_ele] * v[0:n_ele]) + (w[0:n_ele] * w[0:n_ele])) \
       + (inertia * ((ox[0:n_ele] * ox[0:n_ele]) + (oy[0:n_ele] * oy[0:n_ele]) + (oz[0:n_ele] * oz[0:n_ele])))
    ek[0:n_ele] = 0.5 * rho[0:n_ele] * ek[0:n_ele]
    p[0:n_ele] = (Q[0:n_ele,7] - ek) * (k-1.0)
    ####################################
    print >> f_out, "CELL_DATA", n_ele
    print >> f_out, "SCALARS Density float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, rho[i]  
    ####################################       
    print >> f_out, "SCALARS X-Velocity float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, u[i]
    ####################################
    print >> f_out, "SCALARS Y-Velocity float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, v[i]
    ####################################
    print >> f_out, "SCALARS Z-Velocity float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, w[i]
    ####################################       
    print >> f_out, "SCALARS X-Gyration float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, ox[i]
    ####################################
    print >> f_out, "SCALARS Y-Gyration float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, oy[i]
    ####################################
    print >> f_out, "SCALARS Z-Gyration float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, oz[i]
    ####################################
    print >> f_out, "SCALARS Pressure float 1"
    print >> f_out, "LOOKUP_TABLE default"
    for i in range (0, n_ele, 1):
        print >> f_out, p[i]    
    f_out.close()
