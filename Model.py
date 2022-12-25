"""
Date: 20.12.2022
Author: Kaufmann Stefan

verschiedene 2R Robotermodelle  
"""

#from Kinematik_2R import *
from numpy import sin,cos
from Parameter import *
#from Kinematik_2R import f_modell, f_modell_ext


def model_nlin(t,x,u=[0,0],ctr = 'none'):
    """ Nonlinear System Model
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
    #u = eingang(x,t)
    
    q1  = x[0]
    qd1 = x[1]
    q2  = x[2]
    qd2 = x[3]


    if ctr == 'multivariable':
        u = ctr_multi(x)
        

    u1 = u[0]
    u2 = u[1]   


    #qdd = f_modell(qd1,qd2,q1,q2,u1,u2)


    dx[0] = x[1]
    dx[1] = ((I2 + l_s2**2*m2)*(-g*l_s1*m1*cos(q1) - g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) + 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) + u1) + (I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2))/(I1*I2 + I1*l_s2**2*m2 + I2*l1**2*m2 + I2*l_s1**2*m1 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2)
   # dx[1] = qdd[0]    
    dx[2] = x[3]    
    dx[3] = -((I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(-g*l_s1*m1*cos(q1) - g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) + 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2) + u1) + (g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2) - u2)*(I1 + I2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2))/(I1*I2 + I1*l_s2**2*m2 + I2*l1**2*m2 + I2*l_s1**2*m1 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2)
    #dx[3] = qdd[1]



    return dx

    
def model_nlin_ext(t,x,u=[0,0],ctr = 'none'):


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
    #u = eingang(x,t)
    
    q1  = x[0]
    qd1 = x[1]
    q2  = x[2]
    qd2 = x[3]



    #qdd = f_modell_ext(qd1,qd2,q1,q2,u1,u2)

    if ctr == 'multivariable':
        u = ctr_multi_ext(x)
        


    u1 = u[0]
    u2 = u[1]
    #print(u1)

    dx[0] = x[1]
    dx[1] = (R1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(1.0*R2*l_s2*m2*(g*cos(q1 + q2) + l1*qd1**2*sin(q2)) - km2*r2*u2 + qd2*(B2*R2*r2**2 + kb2*km2)) - R2*(I2 + J2*r2**2 + l_s2**2*m2)*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)) - km1*r1*u1 + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1)))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    #dx[1] = qdd[0]
    dx[2] = x[3]
    dx[3] = (R1*(-1.0*R2*l_s2*m2*(g*cos(q1 + q2) + l1*qd1**2*sin(q2)) + km2*r2*u2 - qd2*(B2*R2*r2**2 + kb2*km2))*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + R2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2)*(R1*(g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)) - km1*r1*u1 + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1)))/(R1*R2*(I1*I2 + I1*J2*r2**2 + I1*l_s2**2*m2 + I2*J1*r1**2 + I2*J2*r2**2 + I2*l1**2*m2 + I2*l_s1**2*m1 + J1*J2*r1**2*r2**2 + J1*l_s2**2*m2*r1**2 + J2*l1**2*m2*r2**2 + 2*J2*l1*l_s2*m2*r2**2*cos(q2) + J2*l_s1**2*m1*r2**2 + J2*l_s2**2*m2*r2**2 + l1**2*l_s2**2*m2**2*sin(q2)**2 + l_s1**2*l_s2**2*m1*m2))
    #dx[3] = qdd[1]



    return dx

def ctr_multi_ext(x,u=[0,0]):

    
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


    daq =[0,0]  
    q1  = x[0]
    qd1 = x[1]
    q2  = x[2]
    qd2 = x[3]

    aq1 = u[0]
    aq2 = u[1]


    daq[0] = (R1*(aq1*(I1 + I2 + J1*r1**2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + aq2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)) + qd1*(R1*(B1*r1**2 - 1.0*l1*l_s2*m2*qd2*sin(q2)) + kb1*km1))/R1
    daq[1] = (R2*(aq1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + aq2*(I2 + J2*r2**2 + l_s2**2*m2) + g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2)) + qd2*(B2*R2*r2**2 + kb2*km2))/R2


    return daq


def ctr_multi(x,u=[0,0]):

    
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


    daq =[0,0]  
    q1  = x[0]
    qd1 = x[1]
    q2  = x[2]
    qd2 = x[3]

    aq1 = u[0]
    aq2 = u[1]


    daq[0] = aq1*(I1 + I2 + l1**2*m2 + 2*l1*l_s2*m2*cos(q2) + l_s1**2*m1 + l_s2**2*m2) + aq2*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + g*l_s1*m1*cos(q1) + g*m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) - 1.0*l1*l_s2*m2*qd1*qd2*sin(q2) - 1.0*l1*l_s2*m2*qd2*(qd1 + qd2)*sin(q2)
    daq[1] =  aq1*(I2 + l1*l_s2*m2*cos(q2) + l_s2**2*m2) + aq2*(I2 + l_s2**2*m2) + g*l_s2*m2*cos(q1 + q2) + 1.0*l1*l_s2*m2*qd1**2*sin(q2)


    return daq