"""
Date: 20.12.2022
Author: Kaufmann Stefan

Parameter des Roboters
"""
import numpy as np

# Parameter des Roboters
l1 =   1
l2 =   1
l3 =   1
I1 =   10
I2 =   10
I3 =   10
m1 =   50
m2 =   50
m3 =   50
l_s1 = 0.5
l_s2 = 0.5
l_s3 = 0.5
g =    9.81
J1 =   0.01
J2 =   0.01
J3 =   0.01

# Erweitertes Modell
B1 = 0.01
B2 = 0.01
B3 = 0.01
R1 = 10
R2 = 10
R3 = 10
km1 = 20
km2 = 20
km3 = 20
kb1 = 20
kb2 = 20
kb3 = 20
r1 = 100
r2 = 100
r3 = 100



# define the discretization points
t_start = 0
t_stop = 10
dt = 1e-3

t_sim=np.linspace(t_start, t_stop, int((t_stop - t_start) / dt + 1))
