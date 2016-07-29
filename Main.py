import time
start_time = time.time()
import Input as IN
import numpy as np
import Jacobian as jp
import ReadIn as RI
import ElementBasicInfo as EBI
import CellValue as CV
import Flux
import WallValue as WV
import Output as OP
import math

# Read in input parameter
f_a, f_b, f_c, rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf,\
inertia, dt, TotalStep, RK_Mode, Vis_Mode, Mu, Kappa, Alpha_T, Alpha, Beta, Gamma, TherCon, k, re_mode, Re = IN.input_parameter()

# (1) numper of points; (2) Coordinates of pointes; (3) number of boundary conditions; 
# (4) number of volume elements; (5) boundary conditions; (6) connectivity for volume elements
n_pt, coord = RI.coordinate(f_a)
n_bound, n_ele, bound, ijk = RI.boundary_connectivity(f_b, f_c)

# Calculate (1)volume, (2)center of the volume, (3)areas of all surfaces, (4)center coordinate of all surafces 
#       and (5)all normal vectors of all surfaces of all colume elements
Ele_Vol, Ele_Coor, \
Ele_E_Area, Ele_E_Coor, Ele_E_Nor, \
Ele_W_Area, Ele_W_Coor, Ele_W_Nor, \
Ele_F_Area, Ele_F_Coor, Ele_F_Nor, \
Ele_B_Area, Ele_B_Coor, Ele_B_Nor, \
Ele_N_Area, Ele_N_Coor, Ele_N_Nor, \
Ele_S_Area, Ele_S_Coor, Ele_S_Nor = EBI.element_info(coord, n_ele, ijk)

# Neighbor list or boundary surfaces for all volumes
Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei = EBI.neighbor(n_ele, ijk)
Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei = EBI.boundary_condition(n_ele, n_bound, ijk, bound, \
                                                                                          Ele_E_Nei, Ele_W_Nei, \
                                                                                          Ele_F_Nei, Ele_B_Nei, \
                                                                                          Ele_N_Nei, Ele_S_Nei)

if re_mode == 0:                                                                                          
# Initialize field variables, Q, for all volume elements
    Q = CV.Initiate (rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, n_ele, k)
elif re_mode == 1:
# Restart from the previous simulation
    Q = CV.Restart(n_ele, Re)
else:
    print "Error. No initiation nor restart file"

# Print out the first frame before time marches on
i_step = 1
OP.print_output(Q, n_ele, inertia, k, i_step, coord, n_pt, ijk)



# Time Marching Runge-Kutta Coefficient
stage, RK_a, RK_b, RK_c = EBI.TimeMarching(RK_Mode)
RK_k = np.zeros(shape=(stage, n_ele, 8), dtype=float)

# Time Marching on Here
for i_step in range(2,TotalStep+1,1):
# Runge-Kutta Marching Algorithm Begins
    for j in range(0,stage,1):
        for m in range(0,j+1,1):
# Compute inviscid fluc for all volume elements
                F_inv, G_inv, H_inv, = Flux.compute_inv_flux(n_ele, Q + (dt * RK_a[j, k]*RK_k[m,0:n_ele, 0:8]), inertia, k)

                RK_k[j, 0:n_ele, 0:8] = Flux.Total_Inv_Flux(Ele_Coor, Ele_E_Nei, Ele_W_Nei, Ele_F_Nei, Ele_B_Nei, Ele_N_Nei, Ele_S_Nei, F_inv, G_inv, H_inv, \
                                        Ele_E_Nor, Ele_W_Nor, Ele_F_Nor, Ele_B_Nor, Ele_N_Nor, Ele_S_Nor,\
                                        Ele_E_Area,Ele_W_Area,Ele_F_Area,Ele_B_Area,Ele_N_Area,Ele_S_Area, Ele_Vol, n_ele,\
                                        Ele_E_Coor,Ele_W_Coor,Ele_F_Coor,Ele_B_Coor,Ele_N_Coor,Ele_S_Coor,\
                                        rho, p_inf, u_inf, v_inf, w_inf, ox_inf, oy_inf, oz_inf, inertia, k, Q + (dt * RK_a[j, k]*RK_k[m,0:n_ele, 0:8]))

#                print "There is no much mode for viscous (Vis_Mode = 1) nor invisicid (Vis_Mode = 0). The current value of Vis_Mode is", Vis_Mode


    for l in range(0,stage,1):
        Q[0:n_ele, 0:8] = Q[0:n_ele, 0:8] + (dt * RK_b[l] * RK_k[l, 0:n_ele, 0:8])
# Runge-Kutta Marching Algorithm Ends

# How often does the restart file get generated. Default is every 1000 steps.
    if math.fmod(i_step-1,1000) == 0:
        OP.print_restart(Q, n_ele, Re)
        OP.print_output(Q, n_ele, inertia, k, i_step, coord, n_pt, ijk)
# Time Indicator
    print "===================================================================="
    print "This is the", i_step, "-th step. It has been", (time.time() - start_time), "seconds"        
 
print("--- %s seconds ---" % (time.time() - start_time))
