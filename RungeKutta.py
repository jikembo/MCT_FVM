import numpy as np
import matplotlib.pyplot as plt
import math

def Heun ():
    """
    Heun's Method for time marching. It is a second order two stage explicit method
    repeated index means Einstein summation
    dy/dt = f(y, T)
    y(T + DeltaT) = y(T) + DeltaT * b_i k_i
    k_i = f(T + c_i * DeltaT, y(T) + DeltaT * a_ij k_j)
    """
    stage = 2
    H_a = np.ndarray(shape = (stage, stage), dtype = float)
    H_a[0,0] = 0
    H_a[0,1] = 0
    H_a[1,0] = 1.0
    H_a[1,1] = 0
    
    H_b = np.zeros(stage, dtype = float)
    H_b[0] = 0.5
    H_b[1] = 0.5
    
    H_c = np.zeros(stage, dtype = float)
    H_c[0] = 0
    H_c[1] = 1.0
    return stage, H_a, H_b, H_c




def Kutta():
    """
    Kutta's Method for time marching. It is a third order three stage explicit method
    repeated index means Einstein summation
    dy/dt = f(y, T)
    y(T + DeltaT) = y(T) + DeltaT * b_i k_i
    k_i = f(T + c_i * DeltaT, y(T) + DeltaT * a_ij k_j)
    """
    stage = 3
    H_a = np.ndarray(shape = (stage, stage), dtype = float)
    H_a[0,0] = 0
    H_a[0,1] = 0
    H_a[0,2] = 0
    H_a[1,0] = 0.5
    H_a[1,1] = 0
    H_a[1,2] = 0
    H_a[2,0] = -1.0
    H_a[2,1] = 2.0
    H_a[2,2] = 0
    
    H_b = np.zeros(stage, dtype = float)
    H_b[0] = 1.0/6.0
    H_b[1] = 2.0/3.0
    H_b[2] = 1.0/6.0
    
    H_c = np.zeros(stage, dtype = float)
    H_c[0] = 0
    H_c[1] = 0.5
    H_c[2] = 1.0
    return stage, H_a, H_b, H_c
    
    
    
def Simpson():
    """
    3/8-rule fourth-order Method for time marching. It is a fourth order four stage explicit method
    repeated index means Einstein summation
    dy/dt = f(y, T)
    y(T + DeltaT) = y(T) + DeltaT * b_i k_i
    k_i = f(T + c_i * DeltaT, y(T) + DeltaT * a_ij k_j)
    """
    stage = 4
    H_a = np.ndarray(shape = (stage, stage), dtype = float)
    H_a[0,0] = 0.0
    H_a[0,1] = 0.0
    H_a[0,2] = 0.0
    H_a[0,3] = 0.0
    H_a[1,0] = 1.0/3.0
    H_a[1,1] = 0.0
    H_a[1,2] = 0.0
    H_a[1,3] = 0.0
    H_a[2,0] = -1.0/3.0
    H_a[2,1] = 1.0
    H_a[2,2] = 0.0
    H_a[2,3] = 0.0
    H_a[3,0] = 1.0
    H_a[3,1] = -1.0
    H_a[3,2] = 1.0
    H_a[3,3] = 0.0
    
    H_b = np.zeros(stage, dtype = float)
    H_b[0] = 0.125
    H_b[1] = 0.375
    H_b[2] = 0.375
    H_b[3] = 0.125
    
    H_c = np.zeros(stage, dtype = float)
    H_c[0] = 0.0
    H_c[1] = 1.0/3.0
    H_c[2] = 2.0/3.0
    H_c[3] = 1.0
    return stage, H_a, H_b, H_c  
    
    
      
    
# Runge Kutta Test Session
if __name__ == "__main__":
# Solving dy/dt = -y as an example

    b = np.zeros(3, dtype = float)
    a = np.ndarray(shape = (3,3), dtype = float)
    n, a, b, c = Heun()

    #n, a, b, c = Kutta()

    #n, a, b, c = Simpson()

#    print n, a, b, c

    k = np.zeros(n, dtype = float)

    deltaT = 0.25
    TotalT = 10
    Answer = np.zeros(40, dtype = float)
    Answer[0] = 1

    for i in range(1, 40, 1):
        for j in range(0, n, 1):
            y = 0
            for m in range(0, n, 1):
                y = y + (deltaT * a[j][m] * k[m])
            y = y + Answer[i-1]
            k[j] = y     # The actual function...
        
        for l in range(0, n, 1):
            Answer[i] = Answer[i] + (deltaT * b[l] * k[l])
        Answer[i] = Answer[i-1] + Answer[i]
    
    t = np.arange(0., TotalT, deltaT)      
    plt.plot(t, Answer, 'r--', t, np.exp(t), 'o')
    plt.show()
