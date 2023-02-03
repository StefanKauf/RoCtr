"""
Date: 20.12.2022
Author: Kaufmann Stefan

Parameter des Roboters
"""

######## ***************************************  
## einfaches Modell
######## *************************************** 

# Armlängen  [m]
l1 =   1
l2 =   1
l3 =   0.5
# Abstand zum Massenmittelpunkt vom Gelnk aus [m]
l_s1 = l1/2
l_s2 = l2/2
l_s3 = l3/2  
# Trägheitsmomente [kgm²]
I1 =   10
I2 =   10
I3 =   5      
# Massen  [kg]
m1 =   50
m2 =   50
m3 =   25    
# Erdbeschleunigung [m/s²]
g =    9.81   


######## ***************************************  
## Erweitertes Modell
######## *************************************** 
# Trägheitskonstanten
J1 =   0.01
J2 =   0.01
J3 =   0.01

# Reibkoeffizienten  [Nms/rad]
B1 = 0.01
B2 = 0.01
B3 = 0.01
# Ankerwiderstand [Ohm]
R1 = 10
R2 = 10
R3 = 10
# Motorkonstanten  [Nm/A]
km1 = 20
km2 = 20
km3 = 20
# Motorkonstanten [Vs/rad]
kb1 = 20
kb2 = 20
kb3 = 20

# Übersetzung
r1 = 100
r2 = 100
r3 = 100


