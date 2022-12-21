"""
Date: 20.12.2022
Author: Kaufmann Stefan

verschiedene 2R Robotermodelle  
"""

from Kinematik_2R import *





def model_nlin(x,t,u=[0,0]):
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

    u1 = u[0]
    u2 = u[1]

    qdd = f_modell(qd1,qd2,q1,q2,u1,u2)


    dx[0] = x[1]
    dx[1] = qdd[0]
    dx[2] = x[3]
    dx[3] = qdd[1]



    return dx


    
def model_nlin_ext(x,t,u=[0,0]):
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

    u1 = u[0]
    u2 = u[1]

    qdd = f_modell_ext(qd1,qd2,q1,q2,u1,u2)


    dx[0] = x[1]
    dx[1] = qdd[0]
    dx[2] = x[3]
    dx[3] = qdd[1]



    return dx