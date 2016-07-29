import numpy as np
import WallValue as WV

def inv_flux(Qv, inertia, k):
    """
    Calculate inviscid flux
    Return inviscid flux for i-th volume elements
    """
    F_inv = np.zeros(8, dtype=float)
    G_inv = np.zeros(8, dtype=float)
    H_inv = np.zeros(8, dtype=float)
    rho = Qv[0]
    u = Qv[1] / Qv[0]
    v = Qv[2] / Qv[0]
    w = Qv[3] / Qv[0]
    ox = Qv[4] / Qv[0] / inertia
    oy = Qv[5] / Qv[0] / inertia
    oz = Qv[6] / Qv[0] / inertia
    ek = ((u * u) + (v * v) + (w * w)) \
       + (inertia * ((ox * ox) + (oy * oy) + (oz * oz)))
    ek = 0.5 * rho * ek
    p = (Qv[7] - ek) * (k-1.0)
    
    F_inv[0] = - rho * u
    F_inv[1] = - ( p + (rho * u * u))
    F_inv[2] = - rho * u * v
    F_inv[3] = - rho * u * w
    F_inv[4] = - rho * inertia * ox * u
    F_inv[5] = - rho * inertia * oy * u
    F_inv[6] = - rho * inertia * oz * u
    F_inv[7] = - ((rho * Qv[7]) + p) * u
    
    G_inv[0] = - rho * v
    G_inv[1] = - rho * u * v
    G_inv[2] = - ( p + (rho * v * v))
    G_inv[3] = - rho * v * w
    G_inv[4] = - rho * inertia * ox * v
    G_inv[5] = - rho * inertia * oy * v
    G_inv[6] = - rho * inertia * oy * v
    G_inv[7] = - ((rho * Qv[7]) + p) * v
    
    H_inv[0] = - rho * w
    H_inv[1] = - rho * u * w
    H_inv[2] = - rho * v * w
    H_inv[3] = - ( p + (rho * w * w))
    H_inv[4] = - rho * inertia * ox * w
    H_inv[5] = - rho * inertia * oy * w
    H_inv[6] = - rho * inertia * oy * w
    H_inv[7] = - ((rho * Qv[7]) + p) * w
    return F_inv, G_inv, H_inv

##########################################################################################################

def Total_Inv_Flux(Ele_Coor, Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei, F_inv, G_inv, H_inv, \
                         Ele_E_Nor, Ele_W_Nor, Ele_F_Nor, Ele_B_Nor, Ele_N_Nor, Ele_S_Nor,\
                         Ele_E_Area,Ele_W_Area,Ele_F_Area,Ele_B_Area,Ele_N_Area,Ele_S_Area, Ele_Vol, n_ele,\
                         Ele_E_Coor,Ele_W_Coor,Ele_F_Coor,Ele_B_Coor,Ele_N_Coor,Ele_S_Coor,\
                         rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k, Q):
    """
    This function computes the total inviscid flux in a cell
    """
    total = np.zeros(shape=(n_ele,8), dtype=float)
    V_F = np.zeros(8, dtype=float)
    V_G = np.zeros(8, dtype=float)
    V_H = np.zeros(8, dtype=float)
    for i in range(0, n_ele, 1):
        if Ele_E_Nei[i] < 0:
#            V_F, V_G, V_H = WV.InflowBC_inv(rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k)
            V_F, V_G, V_H = WV.HomoBC_inv(Q[i,0:8], inertia, k, rho)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_E_Nei[i], 0:3], Ele_E_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_E_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_E_Nei[i], 0:3], Ele_E_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_E_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_E_Nei[i], 0:3], Ele_E_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_E_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_E_Nor[i, 0]) + (V_G[0:8] * Ele_E_Nor[i, 1])+ (V_H[0:8] * Ele_E_Nor[i, 2])) * Ele_E_Area[i])


        if Ele_W_Nei[i] < 0:     
#            V_F, V_G, V_H = WV.InflowBC_inv(rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k)
            V_F, V_G, V_H = WV.HomoBC_inv(Q[i,0:8], inertia, k, rho)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_W_Nei[i], 0:3], Ele_W_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_W_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_W_Nei[i], 0:3], Ele_W_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_W_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_W_Nei[i], 0:3], Ele_W_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_W_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_W_Nor[i, 0]) + (V_G[0:8] * Ele_W_Nor[i, 1])+ (V_H[0:8] * Ele_W_Nor[i, 2])) * Ele_W_Area[i])


        if Ele_F_Nei[i] < 0:
            V_F, V_G, V_H = WV.HomoBC_inv(Q[i, 0:8], inertia, k, rho)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_F_Nei[i], 0:3], Ele_F_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_F_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_F_Nei[i], 0:3], Ele_F_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_F_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_F_Nei[i], 0:3], Ele_F_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_F_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_F_Nor[i, 0]) + (V_G[0:8] * Ele_F_Nor[i, 1])+ (V_H[0:8] * Ele_F_Nor[i, 2])) * Ele_F_Area[i])

        if Ele_B_Nei[i] < 0:
            V_F, V_G, V_H = WV.HomoBC_inv(Q[i, 0:8], inertia, k, rho)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_B_Nei[i], 0:3], Ele_B_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_B_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_B_Nei[i], 0:3], Ele_B_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_B_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_B_Nei[i], 0:3], Ele_B_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_B_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_B_Nor[i, 0]) + (V_G[0:8] * Ele_B_Nor[i, 1])+ (V_H[0:8] * Ele_B_Nor[i, 2])) * Ele_B_Area[i])

        if Ele_N_Nei[i] < 0:
#            V_F, V_G, V_H = WV.HomoBC_inv(Q[i, 0:8], inertia, k, rho)
            V_F, V_G, V_H = WV.InflowBC_inv(rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_N_Nei[i], 0:3], Ele_N_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_N_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_N_Nei[i], 0:3], Ele_N_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_N_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_N_Nei[i], 0:3], Ele_N_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_N_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_N_Nor[i, 0]) + (V_G[0:8] * Ele_N_Nor[i, 1])+ (V_H[0:8] * Ele_N_Nor[i, 2])) * Ele_N_Area[i])

        if Ele_S_Nei[i] < 0:
            V_F, V_G, V_H = WV.HomoBC_inv(Q[i, 0:8], inertia, k, rho)
        else: 
            V_F[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_S_Nei[i], 0:3], Ele_S_Coor[i, 0:3], F_inv[i, 0:8], F_inv[Ele_S_Nei[i],0:8])
            V_G[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_S_Nei[i], 0:3], Ele_S_Coor[i, 0:3], G_inv[i, 0:8], G_inv[Ele_S_Nei[i],0:8])
            V_H[0:8] = WV.Wall_Value (Ele_Coor[i, 0:3], Ele_Coor[Ele_S_Nei[i], 0:3], Ele_S_Coor[i, 0:3], H_inv[i, 0:8], H_inv[Ele_S_Nei[i],0:8])
        total[i, 0:8] = total[i, 0:8] + (((V_F[0:8] * Ele_S_Nor[i, 0]) + (V_G[0:8] * Ele_S_Nor[i, 1])+ (V_H[0:8] * Ele_S_Nor[i, 2])) * Ele_S_Area[i])
    
        total[i,0:8] = total[i, 0:8] / Ele_Vol[i]
###########################################################
###########################################################
#        print "========================================"
#        print i,"-th element"
#        print total[i, 0:8]
###########################################################
###########################################################
    return total

##########################################################################################################

def compute_inv_flux(n_ele, Qv, inertia, k):
    F_inv = np.ndarray(shape=(n_ele,8), dtype=float)
    G_inv = np.ndarray(shape=(n_ele,8), dtype=float)
    H_inv = np.ndarray(shape=(n_ele,8), dtype=float)
    for i in range(0, n_ele, 1):
        F_inv[i, 0:8], G_inv[i, 0:8], H_inv[i, 0:8] = inv_flux(Qv[i, 0:8], inertia, k)
##########################################################
##########################################################
#        print i+1, "=th volume."
#        print "F_inv", F_inv[i,0:8]
#        print "G_inv", G_inv[i,0:8]
#        print "H_inv", H_inv[i,0:8]
##########################################################
##########################################################
    return F_inv, G_inv, H_inv

##########################################################################################################
    
# Viscous flux not complete yet Dec. 16, 2015
def vis_flux(Qv, inertia, k, i, coord, Ele_E_Nei, Ele_W_Nei, Ele_F_Nei,\
             Ele_B_Nei, Ele_N_Nei, Ele_S_Nei, E_Coor, W_Coor, F_Coor,\
             B_Coor, N_Coor, S_Coor):
    """
    Calculate viscous flux
    Given: Qv: Field variables
           inertia: microinertia
           k: cp / cv
           i: i-volume element
           coord: everyone's coordinate
           Ele_?_Nei: i-th element's neighbor list
           ?_Coor: coordinate of the surface coordinate of i-th volume element
    Return viscous flux for i-th volume elements
    """
    F_vis = np.zeros(8, dtype=float)
    G_vis = np.zeros(8, dtype=float)
    H_vis = np.zeros(8, dtype=float)
    rho = Qv[0]
    u = Qv[1] / Qv[0]
    v = Qv[2] / Qv[0]
    w = Qv[3] / Qv[0]
    ox = Qv[4] / Qv[0] / inertia
    oy = Qv[5] / Qv[0] / inertia
    oz = Qv[6] / Qv[0] / inertia
    ek = ((u * u) + (v * v) + (w * w)) \
       + (inertia * ((ox * ox) + (oy * oy) + (oz * oz)))
    ek = 0.5 * rho * ek
    p = (Qv[7] - ek) * (k-1.0)
#    return F_vis, G_vis, H_vis
##########################################################################################################


