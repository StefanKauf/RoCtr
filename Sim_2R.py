
"""
Date: 20.12.2022
Author: Kaufmann Stefan

Robot Control - Simulation 2R Roboter 
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Model import model_nlin, model_nlin_ext
# set the initial conditions
x0=[-np.pi/2,0,0,np.pi/4]

# define the discretization points
t_start = 0
t_stop = 10
dt = 1e-3

t_sim=np.linspace(t_start, t_stop, int((t_stop - t_start) / dt + 1))


#scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False, events=None, vectorized=False, args=None, **options)

#solutionOde=odeint(model_nlin,x0,t_sim)
solOde=solve_ivp(model_nlin_ext,[t_start,t_stop],x0,method = 'RK45',t_eval = t_sim)


plt.plot(solOde.t, solOde.y[0]*180/np.pi, 'b', label='q_1')
plt.plot(solOde.t, solOde.y[2]*180/np.pi, 'g', label='q_2')
plt.legend(loc='best')
plt.xlabel('time')
#plt.ylabel('q1(t), q2(t)')
plt.grid()
#plt.savefig('simulation.png')
plt.show()


import Animation as anim
q1 = solOde.y[0]
q2 = solOde.y[2]
anim.plot(q1,q2,dt)