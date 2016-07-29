import numpy as np
import math

def derivative(X_1, X_2, V_1, V_2):
    """
    Find the derivative on x on the cell center based on central difference method
    Given: X_W: the x-coordinate in the west side neighboring cell
           X_E: the x-coordinate in the east side neighboring cell
           V_W: the value in the west side neighboring cell
           V_E: the value in the east side neighboring cell

    Return: D_X: the x-derivative of the value in the targeting cell
    """
    Deri = (V_1 - V_2 ) / (X_1 - X_2)
    return Deri
    
##########################################################################################################

def Restart(n_ele, Re):
    f_restart = open(Re)
    Q_read = np.genfromtxt(f_restart, dtype = 'float', usecols =(0,1,2,3,4,5,6,7), skip_header = 0)
    read_size = Q_read.shape[0]     
    f_restart.close()
    print "============================================================"
    if read_size == n_ele:
        print "Element Number Matches. The simulation shall proceed."
    else:
        print "Element Number does not match. The mesh contains", n_ele, "elements while the restart file has", read_size
    print "============================================================"
    return Q_read
##########################################################################################################
    
def Initiate(rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, n_ele, k):
    """
    Initiate the field variables, Q, with the property in the infinity
    
    Return: Q: Field varialbes
    """    
    Qv = np.ndarray(shape=(n_ele,8), dtype=float)
    Qv [0:n_ele, 0] = rho
    Qv [0:n_ele, 1] = rho * u_inf
    Qv [0:n_ele, 2] = rho * v_inf
    Qv [0:n_ele, 3] = rho * w_inf
    Qv [0:n_ele, 4] = rho * inertia * ox_inf
    Qv [0:n_ele, 5] = rho * inertia * oy_inf
    Qv [0:n_ele, 6] = rho * inertia * oz_inf
    ek = ((u_inf * u_inf) + (v_inf * v_inf) + (w_inf * w_inf)) \
       + (inertia * ((ox_inf * ox_inf) + (oy_inf * oy_inf) + (oz_inf * oz_inf)))
    ek = 0.5 * rho * ek
    Qv [0:n_ele, 7] = (p_inf/(k-1.0)) + ek
    return Qv
##########################################################################################################    

if __name__ == "__main__":
    
#   Test session for the derivative
# Coordinate for all the neighboring cells and the targeting cell
    PE = np.zeros(3, dtype=float)
    PE = ([1.0, 0.2, 0.0])
    PW = np.zeros(3, dtype=float)
    PW = ([-0.4, 0.3, 0.0])
    PF = np.zeros(3, dtype=float)
    PF = ([0.0, -1.0, 0.0])
    PB = np.zeros(3, dtype=float)
    PB = ([0.0, 1.0, 0.0])
    PS = np.zeros(3, dtype=float)
    PS = ([0.0, 0.0, -1.0])
    PN = np.zeros(3, dtype=float)
    PN = ([0.0, 0.0, 1.0])
    PC = np.zeros(3, dtype=float)
    PC = ([0.1, 0.25, 0.])

# Value in all neighboring cells and targeting cell
    QE = 1.0
    QW = 3.0
    QF = 1.0
    QB = 3.0
    QS = 1.0
    QN = 3.0
    QC = 5.0

# Calculate the gradient in targeting cell
    NablaQ = np.zeros(3, dtype=float)
    NablaQ[0] = derivative(PE[0], PW[0], QE, QW)
    NablaQ[1] = derivative(PF[1], PB[1], QF, QB)
    NablaQ[2] = derivative(PS[2], PN[2], QS, QN)
    print "all interior sides, the derivative in the center of the targeting cell", NablaQ

# if the targeting cell is the at the boundary, use the targeting coordinate to replace the boundary
    NablaQ = np.zeros(3, dtype=float)
# if the ease side is boundary
    NablaQ[0] = derivative(PC[0], PW[0], QC, QW)
    NablaQ[1] = derivative(PF[1], PB[1], QF, QB)
    NablaQ[2] = derivative(PS[2], PN[2], QS, QN)
    print "east side is a boundary, the derivative in the center of the targeting cell", NablaQ
# if the front side is boundary
    NablaQ[0] = derivative(PE[0], PW[0], QE, QW)
    NablaQ[1] = derivative(PC[1], PB[1], QC, QB)
    NablaQ[2] = derivative(PS[2], PN[2], QS, QN)
    print "front side is a boundary, the derivative in the center of the targeting cell", NablaQ
# if the south side is boundary
    NablaQ[0] = derivative(PE[0], PW[0], QE, QW)
    NablaQ[1] = derivative(PF[1], PB[1], QF, QB)
    NablaQ[2] = derivative(PC[2], PN[2], QC, QN)
    print "south side is a boundary, the derivative in the center of the targeting cell", NablaQ
