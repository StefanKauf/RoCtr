"""
Date: 20.12.2022
Author: Kaufmann Stefan

verschiedene 2R Robotermodelle  
"""

#from Kinematik_2R import *
import sympy as sym
from numpy import sin,cos
from Parameter import *


######## ***************************************  
##         MODELLE    
##  
######## ***************************************  

def model_nlin(t,x,controller=""):
    """ Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,q2,qd1,qd2]
        t:             time as int
        controller:    u_soll... Solltrajektorie  [q1_soll, q2_soll, qd1_soll, qd2_soll, qdd1_soll, qdd2_soll]
                       ctr ... name of the controller
                       k0 .... Gain matrix
                       k1 .... Gain matrix verlocity


        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """
     
    dx =[0,0,0,0]  
    #u = eingang(x,t)
    
    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]

         
    if controller.ctr == 'multivariable':
        u = ctr_multi(t,x,controller)
    else:
        u =[0,0]
        

    u1 = u[0]
    u2 = u[1]   

    dx[0] = qd1
    dx[1] = qd2
    dx[2] = ((I2 + l_s2**2*m2)*(-g*l_s1*m1*cos(q1) - g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) + 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) + u1) + (I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2))/(I1*I2 + I1*l_s2**2*m2 + I2*l1**2*m2 + I2*l_s1**2*m1 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2)
    dx[3] = -((I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(-g*l_s1*m1*cos(q1) - g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) + 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) + u1) + (g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2)*(I1 + I2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2))/(I1*I2 + I1*l_s2**2*m2 + I2*l1**2*m2 + I2*l_s1**2*m1 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2)
   



    return dx

    
def model_nlin_ext(t,x,controller):


    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        u_soll:        Solltrajektorie  []


        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """
    
    dx =[0,0,0,0]  
  
    
    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]
    

    if controller.trn == True and controller.ctr == 'multivariable':
        u = ctr_multi_ext_transf(t,x,controller,dx)   
        
    elif controller.ctr == 'multivariable' and controller.trn == False:
        u = ctr_multi_ext(t,x,controller)
    else:
        u =[0,0]
        


    u1 = u[0]
    u2 = u[1]
    #print(u1)

    dx[0] = qd1
    dx[1] = qd2   
    dx[2] = -(-R1*(R2*(g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2) + qd2*(B2*R2*r2**2 + kb2*km2))*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + R2*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) - u1) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))*(I2 + J2*r2**2 + l_s2**2*m2))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    dx[3] = (R1*(R2*(-g*l_s2*m2*cos(q1 + q2) - 1.0*l1*l_s2*m2*qd1**2*sin(q2) + u2) - qd2*(B2*R2*r2**2 + kb2*km2))*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + R2*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) - u1) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    
  
    return dx

def model_nlin_(t,x,u):


    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        u_soll:        Solltrajektorie  []


        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """
    
    dx =[0,0,0,0]  
  
    
    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]
    

    u1 = u[0]
    u2 = u[1]
   

    dx[0] = qd1
    dx[1] = qd2   
    dx[2] = -(-R1*(R2*(g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2) + qd2*(B2*R2*r2**2 + kb2*km2))*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + R2*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) - u1) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))*(I2 + J2*r2**2 + l_s2**2*m2))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    dx[3] = (R1*(R2*(-g*l_s2*m2*cos(q1 + q2) - 1.0*l1*l_s2*m2*qd1**2*sin(q2) + u2) - qd2*(B2*R2*r2**2 + kb2*km2))*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + R2*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) - u1) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    
  
    return dx


######## ***************************************  
##         REGLER   
##  für das erweiterte und das nichtlineare Modell
######## ***************************************  
def ctr_multi_ext_transf(t,x,ctr,dx):

    
    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        u_soll:        Solltrajektorie  []


        Returns
        --------
        daq:       input control vektor [daq1, daq2]     
                
    """
    

    u = [0,0]  
    dx = model_nlin_(t,x,u)
    
    u_soll = [0,0,0,0,0,0]

    for i in range(6):
        u_soll[i] = np.interp(t,ctr.t,ctr.ax[i,:])   # Interpolierung von u, damit die Funktion mit ODE-Solver aufgerufen wird kann


    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]

    x_soll =   sym.Matrix(u_soll[0:2])
    xd_soll =  sym.Matrix(u_soll[2:4])
    xdd_soll = sym.Matrix(u_soll[4:6])

    
    # Transfromation vom Gelnk in den Arbeitsraum
    X = tra_J_to_K(x,dx[2:4]) 
    x_   = sym.Matrix(X[0:2])
    xd_  = sym.Matrix(X[2:4])
    xdd_ = sym.Matrix(X[4:6])

    
    # controller im Arbeitsraum
    ax =  -(xdd_*0 -xdd_soll) -ctr.k0*(x_-x_soll) - ctr.k1*(xd_- xd_soll)
    
    aq = tra_K_to_J(x,ax)
    
    
    aq1 = aq[0]
    aq2 = aq[1]

    
    u[0] = (R1*(aq1*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + aq2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))/R1
    u[1] = (R2*(aq1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + aq2*(I2 + J2*r2**2 + l_s2**2*m2) + g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2)) + qd2*(B2*R2*r2**2 + kb2*km2))/R2
    
    

    return u


def ctr_multi_ext(t,x,ctr):

    
    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        u_soll:        Solltrajektorie  []


        Returns
        --------
        daq:       input control vektor [daq1, daq2]     
                
    """


    u = [0,0]  
    u_soll = [0,0,0,0,0,0]

    for i in range(6):
        u_soll[i] = np.interp(t,ctr.t,ctr.u[i,:])   # Interpolierung von u, damit die Funktion mit ODE-Solver aufgerufen wird kann


    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]

    q =        sym.Matrix([q1, q2])
    qd =       sym.Matrix([qd1, qd2])
    q_soll =   sym.Matrix(u_soll[0:2])
    qd_soll =  sym.Matrix(u_soll[2:4])
    qdd_soll = sym.Matrix(u_soll[4:6])


    # controller
    aq = qdd_soll-ctr.k0*(q-q_soll) - ctr.k1*(qd- qd_soll)
     
    aq1 = aq[0]
    aq2 = aq[1]


    u[0] = (R1*(aq1*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + aq2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))/R1
    u[1] = (R2*(aq1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + aq2*(I2 + J2*r2**2 + l_s2**2*m2) + g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2)) + qd2*(B2*R2*r2**2 + kb2*km2))/R2


    return u


def ctr_multi(t,x,ctr):
    
    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        ctr:           u_soll, k0, k1      controller


        Returns
        --------
        daq:       input control vektor [daq1, daq2]     
                
    """
    
    u = [0,0]  

    u_soll = [0,0,0,0,0,0]
    for i in range(6):
        u_soll[i] = np.interp(t,ctr.t,ctr.u[i,:])   # Interpolierung von u, damit die Funktion mit ODE-Solver aufgerufen wird kann


    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]

    q =        sym.Matrix([q1, q2])
    qd =       sym.Matrix([qd1, qd2])
    q_soll =   sym.Matrix(u_soll[0:2])
    qd_soll =  sym.Matrix(u_soll[2:4])
    qdd_soll = sym.Matrix(u_soll[4:6])


    # controller
    aq = qdd_soll-ctr.k0*(q-q_soll) - ctr.k1*(qd- qd_soll)
     
    aq1 = aq[0]
    aq2 = aq[1]

    u[0] = aq1*(I1 + I2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + aq2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)
    u[1] =  aq1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + aq2*(I2 + l_s2**2*m2) + g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2)


    return u



######## ***************************************  
##         TTRANSFORMATION    
##  Kartesicher Raum in den Gelenlswinkelraum
######## ***************************************  

def Ja(x):
    """ Transformation from the Joinspace velocity to kartesian verlocity
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]      
        
        
        Returns
        --------
        (ax:       Kartesian velocitys [x_dot, y_dot, z_dot, omega_X, omega_y, omega_z]  )
        Ja:       analytische Jacobimatrix    
                
    """

    q1  = x[0]
    q2  = x[1]
    '''
    Ja = np.Matrix([[-l1*sin(q1) + l_s2*(-sin(q1)*cos(q2) - sin(q2)*cos(q1)), l_s2*(-sin(q1)*cos(q2) - sin(q2)*cos(q1))],
                    [ l1*cos(q1) + l_s2*(-sin(q1)*sin(q2) + cos(q1)*cos(q2)), l_s2*(-sin(q1)*sin(q2) + cos(q1)*cos(q2))],
                    [                                                      0,                                         0],
                    [                                                      0,                                         0],
                    [                                                      0,                                         0],
                    [                                                      1,                                         1]])

    '''

    Ja = np.matrix([[-l1*sin(q1) - l_s2*sin(q1 + q2), -l_s2*sin(q1 + q2)],
                    [ l1*cos(q1) + l_s2*cos(q1 + q2),  l_s2*cos(q1 + q2)]])
    
       

    return Ja

def Ja_inv(x):
    """ Transformation from the kartesian verlocity to Joinspace velocity 
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]   
        ax:            Kartesian velocitys [x_dot, y_dot, z_dot, omega_X, omega_y, omega_z]    
        
        
        Returns
        --------
        qd:       Join velocitys [q1_dot, q2_dot]     
                
    """

    q1  = x[0]
    q2  = x[1]

    lamda = 0.1*0  # Dämpfung nahe von Singularitäten
    
    Ja_t = sym.Matrix([ [(l1*l_s2**2*(sin(q1 + 2*q2) + sin(3*q1 + 2*q2))/4 - l1*l_s2**2*sin(q1)*cos(q1 + q2)**2 - l1*lamda**2*sin(q1) - l1*sin(q1) - l_s2*lamda**2*sin(q1 + q2))/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2),  (l1**2*l_s2*(sin(q1 - q2) + sin(3*q1 + q2))/4 - l1**2*l_s2*sin(q1 + q2)*cos(q1)**2 - l1*l_s2**2*(sin(q1 + 2*q2) + sin(3*q1 + 2*q2))/4 + l1*l_s2**2*sin(q1)*cos(q1 + q2)**2 + l1*sin(q1) - l_s2*lamda**2*sin(q1 + q2))/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2)],
                        [(-l1*l_s2**2*(cos(q1 + 2*q2) - cos(3*q1 + 2*q2))/4 + l1*l_s2**2*sin(q1 + q2)**2*cos(q1) + l1*lamda**2*cos(q1) + l1*cos(q1) + l_s2*lamda**2*cos(q1 + q2))/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2), (-l1**2*l_s2*(cos(q1 - q2) - cos(3*q1 + q2))/4 + l1**2*l_s2*sin(q1)**2*cos(q1 + q2) + l1*l_s2**2*(cos(q1 + 2*q2) - cos(3*q1 + 2*q2))/4 - l1*l_s2**2*sin(q1 + q2)**2*cos(q1) - l1*cos(q1) + l_s2*lamda**2*cos(q1 + q2))/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2)],
                        [                                                                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                                                                                                                                    0],
                        [                                                                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                                                                                                                                    0],
                        [                                                                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                                                                                                                                    0],
                        [                                                                                                                           (-l1*l_s2*cos(q2) + lamda**2)/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2),                                                                                                                                                                                   (l1**2 + l1*l_s2*cos(q2) + lamda**2)/(l1**2*l_s2**2*sin(q2)**2 + l1**2*lamda**2 + l1**2 + 2*l1*l_s2*lamda**2*cos(q2) + 2*l_s2**2*lamda**2 + lamda**4 + 2*lamda**2)]])
    # landa stellt einen Zusätzlichen Freiheitsgrad dar mit welchem die Anforderungen maximiert werden können [S.42ff, Skript Automatisierungs-und Regelungstechnik WS2021/22 TU-Wien]  --> zB Dämpfung der Geschwindikgeit im Gelenksraum bei singularitäten
    

    
   
    '''
    Ja_inv = np.matrix([[ cos(q1 + q2)/(l1*sin(q2)),                           sin(q1 + q2)/(l1*sin(q2))],
                    [-(l1*cos(q1) + l_s2*cos(q1 + q2))/(l1*l_s2*sin(q2)), -(l1*sin(q1) + l_s2*sin(q1 + q2))/(l1*l_s2*sin(q2))]])

    '''

    Ja_t[0,0] = cos(q1 + q2)/(l1*sin(q2))
    Ja_t[1,0] = -(l1*cos(q1) + l_s2*cos(q1 + q2))/(l1*l_s2*sin(q2))
    Ja_t[0,1] = sin(q1 + q2)/(l1*sin(q2))
    Ja_t[1,1] =  -(l1*sin(q1) + l_s2*sin(q1 + q2))/(l1*l_s2*sin(q2))

    return Ja_t.T

def Ja_diff(x):
    """ Transformation from the Joinspace velocity to kartesian acelaration 
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]   
        qdd:           steady state [qdd1, qdd2 ]
        
        
        Returns
        --------
        xdd:       Join velocitys [x_dot_dort, y_dot_dot]     
                
    """

    q1  = x[0]
    q2  = x[1]
    qd1 = x[2]
    qd2 = x[3]


    Ja_diff = sym.Matrix([[-l1*qd1*cos(q1) - l_s2*qd1*cos(q1 + q2) - l_s2*qd2*cos(q1 + q2), -l_s2*(qd1 + qd2)*cos(q1 + q2)],
                            [-l1*qd1*sin(q1) - l_s2*qd1*sin(q1 + q2) - l_s2*qd2*sin(q1 + q2), -l_s2*(qd1 + qd2)*sin(q1 + q2)],
                            [                                                              0,                              0],
                            [                                                              0,                              0],
                            [                                                              0,                              0],
                            [                                                              0,                              0]])

    

    return Ja_diff



######## ***************************************  
##         TTRANSFORMATION + Regler    
######## ***************************************  

def tra_K_to_J(x,ax):
    """ Inverse 
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        ctr:           u_soll, k0, k1      controller


        Returns
        --------
        daq:       input control vektor [aq1, aq2]     
                
    """

    axx = sym.Matrix.vstack(x[2:4],sym.zeros(4,1))
    
    
    aq = Ja_inv(x)*((ax-np.dot(Ja_diff(x).T,axx)))


    return aq


def tra_J_to_K(x,qdd):
    """ Inverse 
        Params
         --------
        x:             steady states as [q1,qd1,q2,qd2]
        t:             time as int
        qdd:           [qdd1, qdd2]


        Returns
        --------
        qx:          qx = [x,y, xd, yd, xdd, ydd]   
                
    """

    qx = [0,0,0,0,0,0]
    q1 = x[0]
    q2 = x[1]

    qxx = sym.Matrix.vstack(sym.Matrix(x[2:4]),sym.zeros(4,1))

    T02 = np.matrix([[l1*cos(q1) + l_s2*(-sin(q1)*sin(q2) + cos(q1)*cos(q2))],
                    [ l1*sin(q1) + l_s2*(sin(q1)*cos(q2) + sin(q2)*cos(q1))]])

    qx[0] = T02[0,0]                       # Position    
    qx[1] = T02[1,0]
    v     = np.dot(Ja(x),x[2:4])              # Geschwindigkeit
    qx[2] = v[0,0]
    qx[3] = v[0,1]
    a     = np.dot(Ja(x),qdd).T + np.dot(Ja_diff(x).T,qxx)    # Beschleunigung
    qx[4] = a[0,0]
    qx[5] = a[1,0]
    

    return qx

