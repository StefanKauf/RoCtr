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
from IPython.display import HTML


'''
def plot(q1,q2):
    """ Animation
        Params
         --------
        q1:            first joint angle     as vector
        q2:            second joint angle    as vector
        t:             Time vector


        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """


    history_len = 500  # how many trajectory points to display

    L = l1+l2
    x1 = np.cos(q1)*l1
    y1 = np.sin(q2)*l1

    x2 = np.cos(q1+q2)*l2 +x1
    y2 = np.sin(q1+q2)*l2 +y1

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    trace, = ax.plot([], [], '.-', lw=1, ms=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)






    def animate(i):
        thisx = [0, x1[i], x2[i]]
        thisy = [0, y1[i], y2[i]]

        if i == 0:
            history_x.clear()
            history_y.clear()

        history_x.appendleft(thisx[2])
        history_y.appendleft(thisy[2])

        line.set_data(thisx, thisy)
        trace.set_data(history_x, history_y)
        time_text.set_text(time_template % (i*dt))
        return line, trace, time_text


    ani = animation.FuncAnimation(
        fig, animate, len(y1), interval=dt*1000, blit=True)
    plt.show()
'''

    



def plot(q1,q2):
    """ Animation
        Params
         --------
        q1:            first joint angle     as vector
        q2:            second joint angle    as vector
        t:             Time vector


        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """
   

    L = l1+l2
    x1 = np.cos(q1)*l1
    y1 = np.sin(q2)*l1

    x2 = np.cos(q1+q2)*l2 +x1
    y2 = np.sin(q1+q2)*l2 +y1

   
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
    ax.set_aspect('equal')
    ax.grid()
    ax.set_title('Animation')

    # create objects that will change in the animation. These are
    # initially empty, and will be given new values for each frame
    # in the animation.
    line1, = ax.plot([], [], 'b', lw=2)     # ax.plot returns a list of 2D line objects
    line2, = ax.plot([], [], 'r', lw=2)
    
    def drawframe(n):
        line1.set_data([0,x1[n]],[0,y1[n]])
        line2.set_data([x1[n],x2[n]],[y1[n],y2[n]])
        return (line1, line2)

    anim = animation.FuncAnimation(fig, drawframe, frames=100, interval=20, blit=True)
    plt.show()
    HTML(anim.to_html5_video())

   