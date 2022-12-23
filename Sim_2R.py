
"""
Date: 20.12.2022
Author: Kaufmann Stefan

Robot Control - Simulation 2R Roboter 
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Model import model_nlin

def eingang(x,t,k=[1,1]):
    """ System Input
        Params
         --------
        x:             steady states as [x1,x2]
        t:             time as int
        k:             Controler Gain k = [k1,k2]   mit k1,k2 > 0       
                              

        Returns
        --------
        u:              System input   
                
    """
    k1 = k[0]
    k2 = k[1]
    x1 = x[0]
    x2 = x[1]

    #u = -k2*(x2+k1)-x1-x1**2  # Version 1
    u = -k2*(x2+k1*(x1**2)) - x1 -x1**2 -2*k1*x1**2*(x2+k1*x1**2)

    return u 

# set the initial conditions
x0=[np.pi/2,0,0,np.pi/4]

# define the discretization points
t_start = 0
t_stop = 10
dt = 1e-3

t_sim=np.linspace(t_start, t_stop, int((t_stop - t_start) / dt + 1))



solutionOde=odeint(model_nlin,x0,t_sim)

plt.plot(t_sim, solutionOde[:, 0], 'b', label='q_1')
plt.plot(t_sim, solutionOde[:, 2], 'g', label='q_2')
plt.legend(loc='best')
plt.xlabel('time')
plt.ylabel('q1(t), q2(t)')
plt.grid()
#plt.savefig('simulation.png')
plt.show()


import Animation as anim
q1 = solutionOde[:, 0]
q2 = solutionOde[:, 2]
anim.plot(q1,q2)