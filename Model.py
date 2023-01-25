"""
Date: 20.12.2022
Author: Kaufmann Stefan

verschiedene 2R Robotermodelle  
"""

#from Kinematik_2R import *
import sympy as sym
from numpy import sin,cos
from Parameter import *
#from Kinematik_2R import f_modell, f_modell_ext


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


    if controller.ctr == 'multivariable':
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



######## ***************************************  
##         REGLER   
##  für das erweiterte und das nichtlineare Modell
######## ***************************************  
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