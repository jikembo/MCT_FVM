def input_parameter():
#--------------------------FILE NAMES------------------------------
    f_a = './MeshFiles/8k_coord.dat'                #Coordinates Files
    f_b = './MeshFiles/8k_bounds.dat'               #Boundary Conditions Files
    f_c = './MeshFiles/8k_connect.dat'              #Connectivity Files
#--------------------------INFLOW PARAMETER------------------------
    rho = 1.0                                      #Rho_inf
    p_inf = 1.0                                    #Pressure_inf
    u_inf = 0.0                                    #U_inf
    v_inf = 0.0                                    #V_inf
    w_inf = 1.0                                    #W_inf
    o_x_inf = 0.0                                  #Omega_X_inf
    o_y_inf = 0.0                                  #Omega_Y_inf
    o_z_inf = 0.0                                  #Omega_Z_inf
#    T_inf = 300.0                                  #Temperature_inf
    inertia = 2.0e-6                               #MicroInertia
#--------------------------dt for FLOW CALCULATION-----------------
    dt = 5.0e-5                                    #DT
    TotalStep = 20001                              #Total Time Steps
    RK_Mode = 1                                    #Runge-Kutta Scheme; 1:Heun(2nd order); 2:Kutta(3rd order); 3:Simpson(4th order)
#--------------------------VISCOUS RELATED-------------------------
    Vis_Mode = 0                                   #Vismode; 0:Inviscid; 1:Viscous
    Mu = 1.0e-14                                   #Mu
    Kappa = 9.0e-5                                 #Kappa
    Alpha_T = 1.0e-3                               #Alpha_Temprature
    Alpha = 0.0                                    #Alpha_Moment Stress
    Beta = 0.0                                     #Beta_Moment Stress
    Gamma = 1.0e-7                                 #Gamma
    TherCon = 1.0e-3                               #Thermal Conductivity
    k = 1.4                                        #cp / cv
#------------------------RESTART MENU-------------------------------
    restart = 0                                    #Restart Mode (0 = no restart; 1 = use restart file)
    Re = './Results/restart.dat'                   #Restart Files
    return f_a, f_b, f_c, rho, p_inf, u_inf, v_inf, w_inf, o_x_inf, o_y_inf, o_z_inf,\
           inertia, dt, TotalStep, RK_Mode, Vis_Mode, Mu, Kappa, Alpha_T, Alpha, Beta, Gamma, TherCon, k, restart, Re
