"""
Date: 20.12.2022
Author: Kaufmann Stefan

Kinematik 3R Roboter
"""


## Bib
import numpy as np
import sympy as  sym
import Parameter as param



######## ***************************************  
## 1. Transformation 3R Kette
######## ***************************************  

q1,q2,q3,a1,a2,a3,d1,d2,d3,alpha1,alpha2,alpha3,l_s1,l_s2,l_s3,l1,l2,l3= sym.symbols("q1 q2 q3 a1 a2 a3 d1 d2 d3 alpha1 alpha2 alpha3 l_s1 l_s2 l_s3 l1 l2 l3")

## For nice looking equations
sym.init_printing()


q = sym.Matrix([q1,q2,q3])

M1 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q1),-sym.sin(q1),0],[0,sym.sin(q1), sym.cos(q1), 0],[0, 0, 0, 1]])
G1 = sym.Matrix([[1, 0, 0, 0],[a1, 1, 0, 0],[0, 0, sym.cos(alpha1),-sym.sin(alpha1)],[d1, 0, sym.sin(alpha1), sym.cos(alpha1)]])
M2 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q2),-sym.sin(q2),0],[0,sym.sin(q2), sym.cos(q2), 0],[0, 0, 0, 1]])
G2 = sym.Matrix([[1, 0, 0, 0],[a2, 1, 0, 0],[0, 0, sym.cos(alpha2),-sym.sin(alpha2)],[d2, 0, sym.sin(alpha2), sym.cos(alpha2)]])
M3 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q3),-sym.sin(q3),0],[0,sym.sin(q3), sym.cos(q3), 0],[0, 0, 0, 1]])
G3 = sym.Matrix([[1, 0, 0, 0],[a3, 1, 0, 0],[0, 0, sym.cos(alpha3),-sym.sin(alpha3)],[d3, 0, sym.sin(alpha3), sym.cos(alpha3)]])
T03 = M1*G1*M2*G2*M3*G3
T02 = M1*G1*M2*G2
T01 = M1*G1


######## ***************************************  
## 2. Aufstellen der Jacobimatrix
######## ***************************************  

# Geschwindigkeit in x,y Richtung
Jv_01 = T01[1:3,0].jacobian(q)
Jv_02 = T02[1:3,0].jacobian(q)    # 1:3  --> Jacobimatrix bis zur 2ten Dimension, 1 Zeile ist die homogenisierung
Jv_03 = T03[1:3,0].jacobian(q)


# Rotationsgeshwindigkeit

Jw_01 = M1*sym.Matrix([0, 0, 0, 1])+M1[:,0]
Jw_02 = Jw_01.row_join(T01*sym.Matrix([0, 0, 0, 1]))
Jw_03 = Jw_02.row_join(T02*sym.Matrix([0, 0, 0, 1]))
Jw_01 = Jw_01.row_join(sym.Matrix([[0,0,0,0],[0,0,0,0]]).T)
Jw_02 = Jw_02.row_join(T02*sym.Matrix([0, 0, 0, 0]))

# erste Zeile Löschen
Jw_01 = Jw_01[1:4,:]
Jw_02 = Jw_02[1:4,:]
Jw_03 = Jw_03[1:4,:]






# Substitueren mit den Parametern

Jv_1 = Jv_01.subs({a1:l_s1,alpha1:0,d1:0})
Jv_2 = Jv_02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})
Jv_3 = Jv_03.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0})


Jw_1 = Jw_01.subs({a1:l_s1,alpha1:0,d1:0})
Jw_2 = Jw_02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})
Jw_3 = Jw_03.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0})




######## ***************************************  
## 3. Aufstellen der D-Matrix
######## ***************************************  

m1,m2,m3,I1,I2,I3,J1,J2,J3,B1,B2,B3,R1,R2,R3,km1,km2,km3,kb1,kb2,kb3,r1,r2,r3 =sym.symbols("m1 m2 m3 I1 I2 I3 J1 J2 J3 B1 B2 B3 R1 R2 R3 km1 km2 km3 kb1 kb2 kb3 r1 r2 r3")

D = m1*Jv_1.T*Jv_1 + Jw_1.T*I1*Jw_1 + m2*Jv_2.T*Jv_2 + Jw_2.T*I2*Jw_2 + m3*Jv_3.T*Jv_3 + Jw_3.T*I3*Jw_3 
D = sym.simplify(D)

######## ***************************************  
## 4. Aufstellen der C-Matrix
######## ***************************************  

qd1, qd2,qd3, qdd1, qdd2, qdd3 = sym.symbols("qd1 qd2 qd3 qdd1 qdd2 qdd3")
qd = sym.Matrix([qd1,qd2,qd3])
qdd = sym.Matrix([qdd1,qdd2, qdd3])


c = sym.MutableDenseNDimArray(np.zeros((3,)*3))

n = 3
for i in range(n):
    for j in range(n):
        for k in range(n):            
            c[i,j,k] = 1/2*(sym.diff(D[k,j],q[i])+ sym.diff(D[k,i],q[j]) - sym.diff(D[i,j],q[k]) )
           


C = sym.zeros(n,n)

for i in range(n):
    for k in range(n):
        for j in range(n):
            C[k,i] = C[k,i] + c[i,j,k]*qd[j]
        
    
C = sym.simplify(C)

######## ***************************************  
## 5. Aufstellen der g-Vektors  --> potentielle Energie
######## ***************************************  
g = sym.symbols("g")

T01s = T01.subs({a1:l_s1,alpha1:0,d1:0})
T02s = T02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})
T03s = T02.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0})

# Potentielle Energie
P = sym.Matrix([0, 0, g, 0]).T*T01s[:,0]*m1 + sym.Matrix([0, 0, g, 0]).T*T02s[:,0]*m2 + sym.Matrix([0, 0, g, 0]).T*T03s[:,0]*m3

gv = sym.simplify(P.jacobian(q))



######## ***************************************  
## 6. Bewegungsgleichung einfaches Modell
######## ***************************************  
u1,u2,u3 = sym.symbols("u1 u2 u3")

# tau = D*qdd + C*qd + gv
#qdd_Modell = sym.Function('qdd_modell')(qd1,qd2,qd3,q1,q2,q3,u1,u2,u3)

#qdd_modell = D.inv()*(C*qd+gv.T)


######## ***************************************  
## 7. Erweitertes Modell
######## ***************************************  

J = sym.Matrix([[J1*r1**2, 0, 0],[0, J2*r2**2, 0],[0, 0, J3*r3**2]])
B = sym.Matrix([[r1**2*B1, 0, 0],[0, B2*r2**2, 0],[0, 0,  B3*r3**2]])
R = sym.Matrix([[km1*kb1/R1, 0, 0],[0, km2*kb2/R2, 0],[0, 0, km3*kb3/R3]])

tau = sym.Matrix([r1*km1/R1*u1, r2*km2/R2*u2, r3*km3/R3*u3])

# tau = (D+J)qdd + (C+B+R)qd +gv
M = sym.simplify(D+J)
#qdd_ext = M.inv()*((C+B+R)*qd + gv.T)


######## ***************************************  
## 8.  Einsetzen der Parameter
######## ***************************************  



#qdd_modell_subs = sym.simplify(qdd_modell.subs({l1:param.l1,l2:param.l2,I1:param.I1,I2:param.I2,m1:param.m1,m2:param.m2,l_s1:param.l_s1,l_s2:param.l_s2,g:param.g,J1:param.J1,J2:param.J2, B1:param.B1, B2:param.B2,R1:param.R1,R2:param.R2,r1:param.r1,r2:param.r2,km1:param.km1,km2:param.km2,kb1:param.kb1,kb2:param.kb2}))
#f_modell = sym.lambdify([qd1, qd2,q1,q2,u1,u2], qdd_modell_subs)

#qdd_ext_subs = sym.simplify(qdd_ext.subs({l1:param.l1,l2:param.l2,l3:param.l3,I1:param.I1,I2:param.I2,I3:param.I3,m1:param.m1,m2:param.m2,m3:param.m3,l_s1:param.l_s1,l_s2:param.l_s2,l_s3:param.l_s3,g:param.g,J1:param.J1,J2:param.J2,J3:param.J3,B1:param.B1,B2:param.B2,B3:param.B3,R1:param.R1,R2:param.R2,R3:param.R3,r1:param.r1,r2:param.r2,r3:param.r3,km1:param.km1,km2:param.km2,km3:param.km3,kb1:param.kb1,kb2:param.kb2,kb3:param.kb3}))
#f_modell_ext = sym.lambdify([qd1, qd2,q1,q2,u1,u2], qdd_ext_subs)



