import numpy as np
import math

def Wall_Value (P1, P2, Pw, V1, V2):
    """
    compute the value at the wall for the cell locating in P1.
    calculation is based on contant value at the cell wall
    given: P1: the location of the main cell for calculation
           P2: the location of the neighboring cell
           Pw: the location of the wall between the main cell and the neighboring cell
           V1: value in the main cell
           V2: value in the neighboring cell
    return: Vw: value at the wall
    """
    l1w = math.sqrt(math.pow(P1[0]-Pw[0], 2) + math.pow(P1[1]-Pw[1], 2) + math.pow(P1[2]-Pw[2], 2))
    l2w = math.sqrt(math.pow(P2[0]-Pw[0], 2) + math.pow(P2[1]-Pw[1], 2) + math.pow(P2[2]-Pw[2], 2))
#    Vw = ((V2 * l1w) + (V1 * l2w)) / (l1w + l2w)
    Vw = l1w / (l1w + l2w) * (V2 - V1) + V1

    return Vw
    
 ##########################################################################################################   

def HomoBC_inv(Qc, inertia, k, rho):
    F_inv = np.zeros(8, dtype=float)
    G_inv = np.zeros(8, dtype=float)
    H_inv = np.zeros(8, dtype=float)
    
    rho = Qc[0]  
    u = Qc[1] / Qc[0]
    v = Qc[2] / Qc[0]
    w = Qc[3] / Qc[0]
    ox = Qc[4] / Qc[0] / inertia
    oy = Qc[5] / Qc[0] / inertia
    oz = Qc[6] / Qc[0] / inertia
    ek = ((u * u) + (v * v) + (w * w)) \
       + (inertia * ((ox * ox) + (oy * oy) + (oz * oz)))
    ek = 0.5 * rho * ek
    p = (Qc[7] - ek) * (k-1.0)  
    
    F_inv[0] = - rho * u
    F_inv[1] = - ( p + (rho * u * u))
    F_inv[2] = - rho * u * v
    F_inv[3] = - rho * u * w
    F_inv[4] = - rho * inertia * ox * u
    F_inv[5] = - rho * inertia * oy * u
    F_inv[6] = - rho * inertia * oz * u
    F_inv[7] = - (p + Qc[7]) * u
   
    G_inv[0] = - rho * v
    G_inv[1] = - rho * u * v
    G_inv[2] = - ( p + (rho * v * v))
    G_inv[3] = - rho * v * w
    G_inv[4] = - rho * inertia * ox * v
    G_inv[5] = - rho * inertia * oy * v
    G_inv[6] = - rho * inertia * oy * v
    G_inv[7] = - (p + Qc[7]) * v
    
    H_inv[0] = - rho * w
    H_inv[1] = - rho * u * w
    H_inv[2] = - rho * v * w
    H_inv[3] = - ( p + (rho * w * w))
    H_inv[4] = - rho * inertia * ox * w
    H_inv[5] = - rho * inertia * oy * w
    H_inv[6] = - rho * inertia * oy * w
    H_inv[7] = - (p + Qc[7]) * w
    
    return F_inv, G_inv, H_inv   

 ##########################################################################################################   

def InflowBC_inv(rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k):
    F_inv = np.zeros(8, dtype=float)
    G_inv = np.zeros(8, dtype=float)
    H_inv = np.zeros(8, dtype=float)
    u = u_inf
    v = v_inf
    w = w_inf
    ox = ox_inf
    oy = oy_inf
    oz = oz_inf
    ek = ((u * u) + (v * v) + (w * w)) \
       + (inertia * ((ox * ox) + (oy * oy) + (oz * oz)))
    ek = 0.5 * rho * ek
    p = p_inf
    
    F_inv[0] = - rho * u
    F_inv[1] = - ( p + (rho * u * u))
    F_inv[2] = - rho * u * v
    F_inv[3] = - rho * u * w
    F_inv[4] = - rho * inertia * ox * u
    F_inv[5] = - rho * inertia * oy * u
    F_inv[6] = - rho * inertia * oz * u
    F_inv[7] = - (p + ek + (p/(k-1))) * u
   
    G_inv[0] = - rho * v
    G_inv[1] = - rho * u * v
    G_inv[2] = - ( p + (rho * v * v))
    G_inv[3] = - rho * v * w
    G_inv[4] = - rho * inertia * ox * v
    G_inv[5] = - rho * inertia * oy * v
    G_inv[6] = - rho * inertia * oy * v
    G_inv[7] = - (p + ek + (p/(k-1))) * v
    
    H_inv[0] = - rho * w
    H_inv[1] = - rho * u * w
    H_inv[2] = - rho * v * w
    H_inv[3] = - ( p + (rho * w * w))
    H_inv[4] = - rho * inertia * ox * w
    H_inv[5] = - rho * inertia * oy * w
    H_inv[6] = - rho * inertia * oy * w
    H_inv[7] = - (p + ek + (p/(k-1))) * w



    
    return F_inv, G_inv, H_inv
    
########################################################################################################## 

def WallBC_inv(Qv, inertia, k, rho):
    F_inv = np.zeros(8, dtype=float)
    G_inv = np.zeros(8, dtype=float)
    H_inv = np.zeros(8, dtype=float)
    rho = Qv[0]
    u = 0
    v = 0
    w = 0
    ox = 0
    oy = 0
    oz = 0
    ek = ((u * u) + (v * v) + (w * w)) \
       + (inertia * ((ox * ox) + (oy * oy) + (oz * oz)))
    ek = 0.5 * rho * ek
    p = 1.0
    
    F_inv[0] = - rho * u
    F_inv[1] = - ( p + (rho * u * u))
    F_inv[2] = - rho * u * v
    F_inv[3] = - rho * u * w
    F_inv[4] = - rho * inertia * ox * u
    F_inv[5] = - rho * inertia * oy * u
    F_inv[6] = - rho * inertia * oz * u
    F_inv[7] = - (p + ek + (p/(k-1))) * u
    
    G_inv[0] = - rho * v
    G_inv[1] = - rho * u * v
    G_inv[2] = - ( p + (rho * v * v))
    G_inv[3] = - rho * v * w
    G_inv[4] = - rho * inertia * ox * v
    G_inv[5] = - rho * inertia * oy * v
    G_inv[6] = - rho * inertia * oy * v
    G_inv[7] = - (p + ek + (p/(k-1))) * v
    
    H_inv[0] = - rho * w
    H_inv[1] = - rho * u * w
    H_inv[2] = - rho * v * w
    H_inv[3] = - ( p + (rho * w * w))
    H_inv[4] = - rho * inertia * ox * w
    H_inv[5] = - rho * inertia * oy * w
    H_inv[6] = - rho * inertia * oy * w
    H_inv[7] = - (p + ek + (p/(k-1))) * w
    return F_inv, G_inv, H_inv
    
##########################################################################################################
#Test session for Wall Value
if __name__ == "__main__":

# Coordinate for the targeting cell (N1), neighboring cell (N2), and the wall (N2)
    N1 = np.zeros(3, dtype=float)
    N1 = ([0., 0., 0.])
    N2 = np.zeros(3, dtype=float)
    N2 = ([1., 0., 0.])
    Nw = np.zeros(3, dtype=float)
    Nw = ([0.3, 0., 0.])

# Value in the targeting cell (NV1), neighboring (NV2) and Compute the value on the wall (Nwv)
    NV1 = 5.
    NV2 = 0.
    NwV = Wall_Value(N1, N2, Nw, NV1, NV2)
    print NwV
 
