import numpy as np
import math 

def volume(vol_x, vol_y, vol_z):
    """
    Calculate the volume of the box constructed by x, y and z.
    Output vol: volume of the cell
           vol_cel: the coordinates of the cell, vector
    """
    vol = 0
    for i in range (1,3,1):
        r = math.pow(-1.0,i) * 0.577350269189626
        for j in range (1,3,1):
            s = math.pow(-1.0,j) * 0.577350269189626
            for k in range (1,3,1):
                t = math.pow(-1.0,k) * 0.577350269189626

     #   Derivative of Shape Function    #
                N1r=-0.125 * (1-s) * (1+t)
                N1s=-0.125 * (1-r) * (1+t)
                N1t= 0.125 * (1-r) * (1-s)
                N2r= 0.125 * (1-s) * (1+t)
                N2s=-0.125 * (1+r) * (1+t)
                N2t= 0.125 * (1+r) * (1-s)
                N3r= 0.125 * (1+s) * (1+t)
                N3s= 0.125 * (1+r) * (1+t)
                N3t= 0.125 * (1+r) * (1+s)
                N4r=-0.125 * (1+s) * (1+t)
                N4s= 0.125 * (1-r) * (1+t)
                N4t= 0.125 * (1-r) * (1+s)
                N5r=-0.125 * (1-s) * (1-t)
                N5s=-0.125 * (1-r) * (1-t)
                N5t=-0.125 * (1-r) * (1-s)
                N6r= 0.125 * (1-s) * (1-t)
                N6s=-0.125 * (1+r) * (1-t)
                N6t=-0.125 * (1+r) * (1-s)
                N7r= 0.125 * (1+s) * (1-t)
                N7s= 0.125 * (1+r) * (1-t)
                N7t=-0.125 * (1+r) * (1+s)
                N8r=-0.125 * (1+s) * (1-t)
                N8s= 0.125 * (1-r) * (1-t)
                N8t=-0.125 * (1-r) * (1+s)
    
    #   component of jacobian matrix    #
                xr=(N1r * vol_x[0]) + (N2r * vol_x[1]) + (N3r * vol_x[2]) + (N4r * vol_x[3]) + (N5r * vol_x[4]) + (N6r * vol_x[5]) + (N7r * vol_x[6]) + (N8r * vol_x[7])
                xs=(N1s * vol_x[0]) + (N2s * vol_x[1]) + (N3s * vol_x[2]) + (N4s * vol_x[3]) + (N5s * vol_x[4]) + (N6s * vol_x[5]) + (N7s * vol_x[6]) + (N8s * vol_x[7])
                xt=(N1t * vol_x[0]) + (N2t * vol_x[1]) + (N3t * vol_x[2]) + (N4t * vol_x[3]) + (N5t * vol_x[4]) + (N6t * vol_x[5]) + (N7t * vol_x[6]) + (N8t * vol_x[7])
                yr=(N1r * vol_y[0]) + (N2r * vol_y[1]) + (N3r * vol_y[2]) + (N4r * vol_y[3]) + (N5r * vol_y[4]) + (N6r * vol_y[5]) + (N7r * vol_y[6]) + (N8r * vol_y[7])
                ys=(N1s * vol_y[0]) + (N2s * vol_y[1]) + (N3s * vol_y[2]) + (N4s * vol_y[3]) + (N5s * vol_y[4]) + (N6s * vol_y[5]) + (N7s * vol_y[6]) + (N8s * vol_y[7])
                yt=(N1t * vol_y[0]) + (N2t * vol_y[1]) + (N3t * vol_y[2]) + (N4t * vol_y[3]) + (N5t * vol_y[4]) + (N6t * vol_y[5]) + (N7t * vol_y[6]) + (N8t * vol_y[7])
                zr=(N1r * vol_z[0]) + (N2r * vol_z[1]) + (N3r * vol_z[2]) + (N4r * vol_z[3]) + (N5r * vol_z[4]) + (N6r * vol_z[5]) + (N7r * vol_z[6]) + (N8r * vol_z[7])
                zs=(N1s * vol_z[0]) + (N2s * vol_z[1]) + (N3s * vol_z[2]) + (N4s * vol_z[3]) + (N5s * vol_z[4]) + (N6s * vol_z[5]) + (N7s * vol_z[6]) + (N8s * vol_z[7])
                zt=(N1t * vol_z[0]) + (N2t * vol_z[1]) + (N3t * vol_z[2]) + (N4t * vol_z[3]) + (N5t * vol_z[4]) + (N6t * vol_z[5]) + (N7t * vol_z[6]) + (N8t * vol_z[7])

    #   jacobian determinant calculation  #
                deter = (xr * ((ys * zt) - (yt * zs))) - (yr * ((xs * zt) - (xt * zs))) + (zr * ((xs * zt) - (xt * zs)))
                vol = vol + deter
    #   calculate the cell center coordinate #
    vol_cel_x=0
    vol_cel_y=0
    vol_cel_z=0
    for i in range(0,8,1):
        vol_cel_x = vol_cel_x + vol_x[i]                    
        vol_cel_y = vol_cel_y + vol_y[i]
        vol_cel_z = vol_cel_z + vol_z[i]
    vol_cel_x = vol_cel_x * 0.125        
    vol_cel_y = vol_cel_y * 0.125
    vol_cel_z = vol_cel_z * 0.125
    
    vol_cel=np.array([vol_cel_x, vol_cel_y, vol_cel_z])
    
    return vol, vol_cel
    
##########################################################################################################

def area(area_x, area_y, area_z):
    """
    calculate normal vector and the area
    Output nv: normal vector of the surface, vector
           area: area of the surface
           a_cel: the coordinates of the surface, vector
    """
    #    calculate area    with  Bretschneider's formula       
    l=np.zeros(4, dtype=float)
    l[0] = math.pow((area_x[1]-area_x[0]),2) + math.pow((area_y[1]-area_y[0]),2) + math.pow((area_z[1]-area_z[0]),2)
    l[0] = math.sqrt(l[0])
    l[1] = math.pow((area_x[2]-area_x[1]),2) + math.pow((area_y[2]-area_y[1]),2) + math.pow((area_z[2]-area_z[1]),2)
    l[1] = math.sqrt(l[1])
    l[2] = math.pow((area_x[3]-area_x[2]),2) + math.pow((area_y[3]-area_y[2]),2) + math.pow((area_z[3]-area_z[2]),2)
    l[2] = math.sqrt(l[2])
    l[3] = math.pow((area_x[3]-area_x[0]),2) + math.pow((area_y[3]-area_y[0]),2) + math.pow((area_z[3]-area_z[0]),2)
    l[3] = math.sqrt(l[3])
    s = (l[0]+l[1]+l[2]+l[3]) * 0.5

    a = np.array([area_x[1]-area_x[0],area_y[1]-area_y[0],area_z[1]-area_z[0]])
    b = np.array([area_x[1]-area_x[2],area_y[1]-area_y[2],area_z[1]-area_z[2]])
    c = np.array([area_x[3]-area_x[2],area_y[3]-area_y[2],area_z[3]-area_z[2]])
    d = np.array([area_x[3]-area_x[0],area_y[3]-area_y[0],area_z[3]-area_z[0]])
    
    a_d = a * d
    s_a_d = a_d[0] + a_d[1] + a_d[2]
    l_a = math.sqrt(math.pow(a[0],2) + math.pow(a[1],2) + math.pow(a[2],2))
    l_d = math.sqrt(math.pow(d[0],2) + math.pow(d[1],2) + math.pow(d[2],2))
    theta_1 = math.acos(s_a_d/l_a/l_d)
    
    b_c = b * c
    s_b_c = b_c[0] + b_c[1] + b_c[2]
    l_b = math.sqrt(math.pow(b[0],2) + math.pow(b[1],2) + math.pow(b[2],2))
    l_c = math.sqrt(math.pow(c[0],2) + math.pow(c[1],2) + math.pow(c[2],2))
    theta_2 = math.acos(s_b_c/l_b/l_c)    
    
    theta = 0.5 * (theta_1 + theta_2)
    
    area = math.sqrt((s-l[0]) * (s-l[1]) * (s-l[2]) * (s-l[3])) - (l[0] * l[1] * l[2] * l[3] * math.pow(math.cos(theta),2.0))
            
    #   calculate normal vector                    
    v_1=np.zeros(3, dtype=float)
    v_2=np.zeros(3, dtype=float)
    nv=np.zeros(3, dtype=float)
    v_1[0] = area_x[2]-area_x[0]
    v_1[1] = area_y[2]-area_y[0]
    v_1[2] = area_z[2]-area_z[0]
    v_2[0] = area_x[3]-area_x[1]
    v_2[1] = area_y[3]-area_y[1]
    v_2[2] = area_z[3]-area_z[1]
    nv[0] =  (v_1[1] * v_2[2]) - (v_2[1] * v_1[2])
    nv[1] = -(v_1[0] * v_2[2]) + (v_2[0] * v_1[2])
    nv[2] =  (v_1[0] * v_2[1]) - (v_2[0] * v_1[1])
    mag = math.sqrt(math.pow(nv[0],2) + math.pow(nv[1],2) + math.pow(nv[2],2))
    nv[0] = nv[0]/mag
    nv[1] = nv[1]/mag
    nv[2] = nv[2]/mag
     
    #   calculate the cell center coordinate       #
    a_cel_x = 0
    a_cel_y = 0
    a_cel_z = 0
    for i in range(0,4,1):
        a_cel_x = a_cel_x + area_x[i] 
        a_cel_y = a_cel_y + area_y[i] 
        a_cel_z = a_cel_z + area_z[i] 
    a_cel_x = a_cel_x * 0.25
    a_cel_y = a_cel_y * 0.25
    a_cel_z = a_cel_z * 0.25
    
    a_cel=np.array([a_cel_x, a_cel_y, a_cel_z])
    
    return nv, area, a_cel
    
    
    
# Jacobian Test Session
if __name__ == "__main__":   
    x=np.zeros(8, dtype=float)
    y=np.zeros(8, dtype=float)
    z=np.zeros(8, dtype=float)
    x[0]=1.
    y[0]=0.
    z[0]=0.
    x[1]=2.
    y[1]=0.
    z[1]=0.
    x[2]=2.
    y[2]=1.
    z[2]=0.
    x[3]=1.
    y[3]=1.
    z[3]=0.
    x[4]=0.
    y[4]=0.
    z[4]=-1.
    x[5]=1.
    y[5]=0.
    z[5]=-1.
    x[6]=1.
    y[6]=1.
    z[6]=-1.
    x[7]=0.
    y[7]=1.
    z[7]=-1.
    answer_vol, cel_coor = volume(x, y, z)
    print "the cell volume is ", answer_vol, "the cell coordinates are", cel_coor

    a_x=np.array([x[7], x[6], x[5], x[4]])
    a_y=np.array([y[7], y[6], y[5], y[4]])
    a_z=np.array([z[7], z[6], z[5], z[4]])
    nor_v=np.zeros(3, dtype=float)
    nor_v, answer_area, area_coor = area (a_x, a_y, a_z)
    print "the normal vector of the surface", nor_v
    print "the area of the surface is ", answer_area, "locating at", area_coor

