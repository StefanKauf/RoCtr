"""
Date: 21.12.2022
Author: Kaufmann Stefan

Animation 2D plot
"""

# Load Robotermodell
from Parameter import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython import display


def plot(q1,q2,dt):
    """ Animation
        Params
         --------
        q1:            first joint angle     as vector
        q2:            second joint angle    as vector
        dt:            Timestep


        Returns
        --------
          
                
    """
   

    L = l1+l2
    x1 = np.cos(q1)*l1
    y1 = np.sin(q1)*l1
    
    x2 = np.add(np.cos(np.add(q1,q2))*l2,x1)
    y2 = np.add(np.sin(np.add(q1,q2))*l2,y1)

   
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
    ax.set_aspect('equal')
    ax.grid()
    ax.set_title('Animation')

    # create objects that will change in the animation. These are
    # initially empty, and will be given new values for each frame
    # in the animation.
    line1, = ax.plot([], [], 'b', lw=4)     # ax.plot returns a list of 2D line objects
    line2, = ax.plot([], [], 'r', lw=4, ms=2)
    
    
    def drawframe(n):        
        line1.set_data([0,x1[n]],[0,y1[n]])
        line2.set_data([x1[n],x2[n]],[y1[n],y2[n]])
        return (line1,line2)

    anim = animation.FuncAnimation(fig, drawframe, frames=x1.size, interval=dt*1000, blit=True)
    video = anim.to_html5_video()
    html = display.HTML(video)
    display.display(html)
    plt.close()                   # avoid plotting a spare static plot



