import numpy as np
import Jacobian as jp
import RungeKutta as RK

def element_info(coord, n_ele, ijk):
    """
    Calculate basic element information
    return Ele_Vol: volume of elements
           Ele_Coor: center coordinate of volume element
           Ele_?_Area: the area of the ?-th side of the element
           Ele_?_Coor: the center coordinate of the ?-th side of the element
           Ele_?_Nor: the normal vector of the ?-th side of the element
           ? = E, W, F, B, N, S
    """
    Ele_Vol = np.zeros(n_ele, dtype = float)
    Ele_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_E_Area = np.zeros(n_ele, dtype = float)
    Ele_E_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_E_Nor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_W_Area = np.zeros(n_ele, dtype = float)
    Ele_W_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_W_Nor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_F_Area = np.zeros(n_ele, dtype = float)
    Ele_F_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_F_Nor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_B_Area = np.zeros(n_ele, dtype = float)
    Ele_B_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_B_Nor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_N_Area = np.zeros(n_ele, dtype = float)
    Ele_N_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_N_Nor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_S_Area = np.zeros(n_ele, dtype = float)
    Ele_S_Coor = np.ndarray(shape=(n_ele,3), dtype=float)
    Ele_S_Nor = np.ndarray(shape=(n_ele,3), dtype=float)    
    X_Vol = np.zeros(8, dtype = float)
    Y_Vol = np.zeros(8, dtype = float)
    Z_Vol = np.zeros(8, dtype = float)
    # Calculate element volume and area on each surface
    print "============================================================"
    print "=============== Element Volume ============================="
    for i in range(0,n_ele, 1):
        X_Vol = ([coord[ijk[i,0]-1,0], coord[ijk[i,1]-1,0], coord[ijk[i,5]-1,0], \
                  coord[ijk[i,4]-1,0], coord[ijk[i,3]-1,0], coord[ijk[i,2]-1,0], \
                  coord[ijk[i,6]-1,0], coord[ijk[i,7]-1,0]])
        Y_Vol = ([coord[ijk[i,0]-1,1], coord[ijk[i,1]-1,1], coord[ijk[i,5]-1,1], \
                  coord[ijk[i,4]-1,1], coord[ijk[i,3]-1,1], coord[ijk[i,2]-1,1], \
                  coord[ijk[i,6]-1,1], coord[ijk[i,7]-1,1]])
        Z_Vol = ([coord[ijk[i,0]-1,2], coord[ijk[i,1]-1,2], coord[ijk[i,5]-1,2], \
                  coord[ijk[i,4]-1,2], coord[ijk[i,3]-1,2], coord[ijk[i,2]-1,2], \
                  coord[ijk[i,6]-1,2], coord[ijk[i,7]-1,2]])
        Ele_Vol[i], Ele_Coor[i,0:3] = jp.volume(X_Vol, Y_Vol, Z_Vol)

        X_Vol = ([coord[ijk[i,2]-1,0],coord[ijk[i,6]-1,0], coord[ijk[i,5]-1,0], coord[ijk[i,1]-1,0]])
        Y_Vol = ([coord[ijk[i,2]-1,1],coord[ijk[i,6]-1,1], coord[ijk[i,5]-1,1], coord[ijk[i,1]-1,1]])
        Z_Vol = ([coord[ijk[i,2]-1,2],coord[ijk[i,6]-1,2], coord[ijk[i,5]-1,2], coord[ijk[i,1]-1,2]])
        Ele_E_Nor[i,0:3], Ele_E_Area[i], Ele_E_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])
    
        X_Vol = ([coord[ijk[i,7]-1,0],coord[ijk[i,3]-1,0], coord[ijk[i,0]-1,0], coord[ijk[i,4]-1,0]])
        Y_Vol = ([coord[ijk[i,7]-1,1],coord[ijk[i,3]-1,1], coord[ijk[i,0]-1,1], coord[ijk[i,4]-1,1]])
        Z_Vol = ([coord[ijk[i,7]-1,2],coord[ijk[i,3]-1,2], coord[ijk[i,0]-1,2], coord[ijk[i,4]-1,2]])
        Ele_W_Nor[i,0:3], Ele_W_Area[i], Ele_W_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])

        X_Vol = ([coord[ijk[i,3]-1,0],coord[ijk[i,2]-1,0], coord[ijk[i,1]-1,0], coord[ijk[i,0]-1,0]])
        Y_Vol = ([coord[ijk[i,3]-1,1],coord[ijk[i,2]-1,1], coord[ijk[i,1]-1,1], coord[ijk[i,0]-1,1]])
        Z_Vol = ([coord[ijk[i,3]-1,2],coord[ijk[i,2]-1,2], coord[ijk[i,1]-1,2], coord[ijk[i,0]-1,2]])
        Ele_F_Nor[i,0:3], Ele_F_Area[i], Ele_F_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])

        X_Vol = ([coord[ijk[i,6]-1,0],coord[ijk[i,7]-1,0], coord[ijk[i,4]-1,0], coord[ijk[i,5]-1,0]])
        Y_Vol = ([coord[ijk[i,6]-1,1],coord[ijk[i,7]-1,1], coord[ijk[i,4]-1,1], coord[ijk[i,5]-1,1]])
        Z_Vol = ([coord[ijk[i,6]-1,2],coord[ijk[i,7]-1,2], coord[ijk[i,4]-1,2], coord[ijk[i,5]-1,2]])
        Ele_B_Nor[i,0:3], Ele_B_Area[i], Ele_B_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])

        X_Vol = ([coord[ijk[i,0]-1,0],coord[ijk[i,1]-1,0], coord[ijk[i,5]-1,0], coord[ijk[i,4]-1,0]])
        Y_Vol = ([coord[ijk[i,0]-1,1],coord[ijk[i,1]-1,1], coord[ijk[i,5]-1,1], coord[ijk[i,4]-1,1]])
        Z_Vol = ([coord[ijk[i,0]-1,2],coord[ijk[i,1]-1,2], coord[ijk[i,5]-1,2], coord[ijk[i,4]-1,2]])
        Ele_N_Nor[i,0:3], Ele_N_Area[i], Ele_N_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])

        X_Vol = ([coord[ijk[i,7]-1,0],coord[ijk[i,6]-1,0], coord[ijk[i,2]-1,0], coord[ijk[i,3]-1,0]])
        Y_Vol = ([coord[ijk[i,7]-1,1],coord[ijk[i,6]-1,1], coord[ijk[i,2]-1,1], coord[ijk[i,3]-1,1]])
        Z_Vol = ([coord[ijk[i,7]-1,2],coord[ijk[i,6]-1,2], coord[ijk[i,2]-1,2], coord[ijk[i,3]-1,2]])
        Ele_S_Nor[i,0:3], Ele_S_Area[i], Ele_S_Coor[i, 0:3] = jp.area(X_Vol[0:4], Y_Vol[0:4], Z_Vol[0:4])
        print i, "-th element: volume -", Ele_Vol[i]
    return Ele_Vol, Ele_Coor, \
        Ele_E_Area, Ele_E_Coor, Ele_E_Nor, \
        Ele_W_Area, Ele_W_Coor, Ele_W_Nor, \
        Ele_F_Area, Ele_F_Coor, Ele_F_Nor, \
        Ele_B_Area, Ele_B_Coor, Ele_B_Nor, \
        Ele_N_Area, Ele_N_Coor, Ele_N_Nor, \
        Ele_S_Area, Ele_S_Coor, Ele_S_Nor 
        
##########################################################################################################
        
def neighbor(n_ele, ijk):
    """
    Find the neighbors of the elements
    return Ele_?_Nei: the element number of the ? side neighbor of the element
           ? = E, W, F, B, N, S
    """
    Ele_E_Nei = np.zeros(n_ele, dtype = int)
    Ele_W_Nei = np.zeros(n_ele, dtype = int)
    Ele_F_Nei = np.zeros(n_ele, dtype = int)
    Ele_B_Nei = np.zeros(n_ele, dtype = int)
    Ele_N_Nei = np.zeros(n_ele, dtype = int)
    Ele_S_Nei = np.zeros(n_ele, dtype = int)
    for i in range (0, n_ele, 1):
        for j in range(0, n_ele, 1):
            if i == j:
                continue
            elif ([ijk[i,2], ijk[i,6], ijk[i,5], ijk[i,1]]) == ([ijk[j,3], ijk[j,7], ijk[j,4], ijk[j,0]]):
    # i-th east side neighbor is j-th west side neighbor
                Ele_E_Nei[i] = j
                Ele_W_Nei[j] = i

    # i-th front side neighbor is j-th back side neighbor
            elif ([ijk[i,0], ijk[i,1], ijk[i,2], ijk[i,3]]) == ([ijk[j,4], ijk[j,5], ijk[j,6], ijk[j,7]]):
                Ele_F_Nei[i] = j
                Ele_B_Nei[j] = i

    # i-th north side neighbor is j-th south side neighbor
            elif ([ijk[i,0], ijk[i,1], ijk[i,5], ijk[i,4]]) == ([ijk[j,3], ijk[j,2], ijk[j,6], ijk[j,7]]):
                Ele_N_Nei[i] = j
                Ele_S_Nei[j] = i
    return Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei
    
##########################################################################################################
    
def boundary_condition(n_ele, n_bound, ijk, bound, Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei):
    """
    Locate all boundaries
    return updated neighbor list
    """
    Comp_1 = np.zeros(4, dtype = int)
    for i in range (0, n_bound, 1):
        for j in range (0, n_ele,1): 
            if ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 2], ijk[j, 1], ijk[j, 5], ijk[j, 6]]):
                Ele_E_Nei[j] = -1 #bound[i,0]
            elif ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 0], ijk[j, 3], ijk[j, 7], ijk[j, 4]]):
                Ele_W_Nei[j] = -2 #bound[i,0]
            elif ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 2], ijk[j, 1], ijk[j, 0], ijk[j, 3]]):
                Ele_F_Nei[j] = -3 #bound[i,0]
            elif ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 5], ijk[j, 6], ijk[j, 7], ijk[j, 4]]):
                Ele_B_Nei[j] = -4 #bound[i,0]
            elif ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 1], ijk[j, 5], ijk[j, 4], ijk[j, 0]]):
                Ele_N_Nei[j] = -5 #bound[i,0]
            elif ([bound[i,1], bound[i,2], bound[i,3], bound[i,4]]) == ([ijk[j, 2], ijk[j, 6], ijk[j, 7], ijk[j, 3]]):
                Ele_S_Nei[j] = -6 #bound[i,0]
    print "============================================================"
    print "============ Element Neighbor List ========================="
    for i in range (0, n_ele, 1):
        print i, "-th element, East:", Ele_E_Nei[i], "West:", Ele_W_Nei[i],\
                              "Front:", Ele_F_Nei[i], "Back:", Ele_B_Nei[i],\
                              "North:", Ele_N_Nei[i], "South", Ele_S_Nei[i]
    return Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei
##########################################################################################################
def TimeMarching(Mode):
    """
    Decide which explicit Runge-Kutta scheme for time marching
    Return the coefficients and number of the stages
    """
    if Mode == 1:
        s, a, b, c = RK.Heun()
        print "Heun's method is chosen. It is a second order algorithm."
    elif Mode == 2:
        s, a, b, c = RK.Kutta()
        print "Kutta's method is chosen. It is a third order algorithm."
    elif Mode == 3: 
        s, a, b, c = RK.Simpson()
        print "Simpson's 3/8 method is chosen. It is a fourth order algorithm."
    else:
        print "ERROR: No Time Maching Algorithm was chosen." 
    print "=========================================================="
    print "============= Time Marching Algorithm Coefficients ======="
    for i in range(0, s, 1):
        print c[i], a[i,0:s]
    print b
    print "It is a", s, "stage R-K Scheme."
     
    return s, a, b, c
##########################################################################################################
