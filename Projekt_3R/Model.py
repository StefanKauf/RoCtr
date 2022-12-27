"""
Date: 20.12.2022
Author: Kaufmann Stefan

verschiedene 2R Robotermodelle  
"""

#from Kinematik_2R import *
import sympy as sym
from numpy import sin,cos
from Parameter import *

    
def model_nlin_ext(t,x,controller):


    """ Extendend Nonlinear System Model
        Params
         --------
        x:             steady states as [q1,q2,q3,qd1,qd2,qd3]
        t:             time as int
        controller:    u_soll... Solltrajektorie  [q1_soll, q2_soll, q3_soll,qd1_soll, qd2_soll, qd3_soll, qdd1_soll, qdd2_soll, qdd3_soll]  vector n = 9
                       ctr ... name of the controller
                       k0 .... Gain matrix 3*3
                       k1 .... Gain matrix verlocity  3*3


        Returns
        --------
        dx:       chance of the state as a vektor [qd1,qd2,qd3,qdd1,qdd2,qdd3]     
                
    """
    
    dx =[0,0,0,0,0,0]  
  
    
    q1  = x[0]
    q2  = x[1]
    q3  = x[2]
    qd1 = x[3]
    qd2 = x[4]
    qd3 = x[5]


    if controller.ctr == 'multivariable':
        u = ctr_multi_ext(t,x,controller)
    else:
        u = controller.u
        


    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    #print(u1)

    dx[0] = qd1
    dx[1] = qd2
    dx[2] = qd3
    dx[3] = ((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2)*(-g*(l_s1*m1*cos(q1) + m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + m3*(l1*cos(q1) + l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))*(qd1 + qd2 + qd3) - qd1*(B1*r1**2 - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3)) + kb1*km1/R1) - qd2*(-1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))) + u1)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + (-(I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) + (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3))*(-g*(l_s2*m2*cos(q1 + q2) + m3*(l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l2*l_s3*m3*qd3*(qd1 + qd2 + qd3)*sin(q3) - qd1*(1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l2*l_s3*m3*qd3*sin(q3)) - qd2*(B2*r2**2 - 1.0*l2*l_s3*m3*qd3*sin(q3) + kb2*km2/R2) + u2)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + ((I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))*(-g*l_s3*m3*cos(q1 + q2 + q3) - 1.0*l2*l_s3*m3*qd2*(qd1 + qd2)*sin(q3) - 1.0*l_s3*m3*qd1*(l2*qd2*sin(q3) + qd1*(l1*sin(q2 + q3) + l2*sin(q3))) - qd3*(B3*r3**2 + kb3*km3/R3) + u3)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))
    dx[4] = (-(I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) + (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3))*(-g*(l_s1*m1*cos(q1) + m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + m3*(l1*cos(q1) + l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))*(qd1 + qd2 + qd3) - qd1*(B1*r1**2 - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3)) + kb1*km1/R1) - qd2*(-1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))) + u1)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + ((I3 + J3*r3**2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2)*(-g*(l_s2*m2*cos(q1 + q2) + m3*(l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l2*l_s3*m3*qd3*(qd1 + qd2 + qd3)*sin(q3) - qd1*(1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l2*l_s3*m3*qd3*sin(q3)) - qd2*(B2*r2**2 - 1.0*l2*l_s3*m3*qd3*sin(q3) + kb2*km2/R2) + u2)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + (-(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))*(-g*l_s3*m3*cos(q1 + q2 + q3) - 1.0*l2*l_s3*m3*qd2*(qd1 + qd2)*sin(q3) - 1.0*l_s3*m3*qd1*(l2*qd2*sin(q3) + qd1*(l1*sin(q2 + q3) + l2*sin(q3))) - qd3*(B3*r3**2 + kb3*km3/R3) + u3)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))
    dx[5] = ((I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))*(-g*(l_s1*m1*cos(q1) + m2*(l1*cos(q1) + l_s2*cos(q1 + q2)) + m3*(l1*cos(q1) + l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))*(qd1 + qd2 + qd3) - qd1*(B1*r1**2 - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3)) + kb1*km1/R1) - qd2*(-1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l1*qd2*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l_s3*m3*qd3*(l1*sin(q2 + q3) + l2*sin(q3))) + u1)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + (-(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))*(-g*(l_s2*m2*cos(q1 + q2) + m3*(l2*cos(q1 + q2) + l_s3*cos(q1 + q2 + q3))) + 1.0*l2*l_s3*m3*qd3*(qd1 + qd2 + qd3)*sin(q3) - qd1*(1.0*l1*qd1*(l2*m3*sin(q2) + l_s2*m2*sin(q2) + l_s3*m3*sin(q2 + q3)) - 1.0*l2*l_s3*m3*qd3*sin(q3)) - qd2*(B2*r2**2 - 1.0*l2*l_s3*m3*qd3*sin(q3) + kb2*km2/R2) + u2)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)) + ((I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2)*(-g*l_s3*m3*cos(q1 + q2 + q3) - 1.0*l2*l_s3*m3*qd2*(qd1 + qd2)*sin(q3) - 1.0*l_s3*m3*qd1*(l2*qd2*sin(q3) + qd1*(l1*sin(q2 + q3) + l2*sin(q3))) - qd3*(B3*r3**2 + kb3*km3/R3) + u3)/((I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) - (I3 + J3*r3**2 + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3)**2 - (I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I1 + I2 + I3 + J1*r1**2 + l1**2*m2 + l1**2*m3 + 2*l1*l2*m3*cos(q2) + 2*l1*l_s2*m2*cos(q2) + 2*l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s1**2*m1 + l_s2**2*m2 + l_s3**2*m3) + 2*(I3 + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)*(I2 + I3 + l1*l2*m3*cos(q2) + l1*l_s2*m2*cos(q2) + l1*l_s3*m3*cos(q2 + q3) + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3) - (I3 + l1*l_s3*m3*cos(q2 + q3) + l2*l_s3*m3*cos(q3) + l_s3**2*m3)**2*(I2 + I3 + J2*r2**2 + l2**2*m3 + 2*l2*l_s3*m3*cos(q3) + l_s2**2*m2 + l_s3**2*m3))

    return dx

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