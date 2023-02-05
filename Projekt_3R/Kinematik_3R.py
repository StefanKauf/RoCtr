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
#sym.init_printing()


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
Jv_01 = T01[1:4,0].jacobian(q)
Jv_02 = T02[1:4,0].jacobian(q)    # 1:3  --> Jacobimatrix bis zur 2ten Dimension, 1 Zeile ist die homogenisierung
Jv_03 = T03[1:4,0].jacobian(q)


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
Jv_3 = sym.simplify(Jv_03.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0}))


Jw_1 = Jw_01.subs({a1:l_s1,alpha1:0,d1:0})
Jw_2 = Jw_02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})
Jw_3 = sym.simplify(Jw_03.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0}))




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
T03s = T03.subs({a1:l1,a2:l2,a3:l_s3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0})

# Potentielle Energie
P = sym.Matrix([0, 0, g, 0]).T*T01s[:,0]*m1 + sym.Matrix([0, 0, g, 0]).T*T02s[:,0]*m2 + sym.Matrix([0, 0, g, 0]).T*T03s[:,0]*m3

gv = sym.simplify(P.jacobian(q))



######## ***************************************  
## 6. Bewegungsgleichung erweitertes Modell
######## ***************************************  
u1,u2,u3 = sym.symbols("u1 u2 u3")
tau = sym.Matrix([u1,u2,u3])


J = sym.Matrix([[J1*r1**2, 0, 0],[0, J2*r2**2, 0],[0, 0, J3*r3**2]])
B = sym.Matrix([[r1**2*B1, 0, 0],[0, B2*r2**2, 0],[0, 0,  B3*r3**2]])
R = sym.Matrix([[km1*kb1/R1, 0, 0],[0, km2*kb2/R2, 0],[0, 0, km3*kb3/R3]])


M = sym.simplify(D+J)
a11,a12,a13,a21,a22,a23,a31,a32,a33 = sym.symbols("a11 a12 a13 a21 a22 a23 a31 a32 a33")

A = sym.Matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
Ainv = sym.simplify(A.inv())

M_inv = Ainv.subs([(a11,M[0,0]),(a12,M[0,1]),(a13,M[0,2]), (a21,M[1,0]), (a22,M[1,1]), (a23,M[1,2]), (a31,M[2,0]), (a32,M[2,1]), (a33,M[2,2]) ])
#M_inv = sym.simplify(M_inv)
qdd_ext = M_inv*(tau-(C+B+R)*qd - gv.T)




######## ***************************************  
## 7. Kinematik
######## *************************************** 

# Vorwärtskinematik
subs_var = [(l1,param.l1),(l2,param.l2),(l3,param.l3),(I1,param.I1),(I2,param.I2),(I3,param.I3),(m1,param.m1),(m2,param.m2),(m3,param.m3),(l_s1,param.l_s1),(l_s2,param.l_s2),(l_s3,param.l_s3),(g,param.g),(J1,param.J1),(J2,param.J2),(J3,param.J3),(B1,param.B1),(B2,param.B2),(B3,param.B3),(R1,param.R1),(R2,param.R2),(R3,param.R3),(r1,param.r1),(r2,param.r2),(r3,param.r3),(km1,param.km1),(km2,param.km2),(km3,param.km3),(kb1,param.kb1),(kb2,param.kb2),(kb3,param.kb3)]
T03sub = T03.subs({a1:l1,a2:l2,a3:l3,alpha1:0, alpha2:0,alpha3:0,d1:0,d2:0,d3:0})



######## ***************************************  
## 9.  analytische Jacobimatrix
######## ***************************************  
Psi,Theta, lamda = sym.symbols("Psi Theta lamda")
B_a = sym.Matrix([[sym.cos(Psi)*sym.cos(Theta),-sym.sin(Theta), 0],[sym.sin(Psi)*sym.sin(Theta), sym.cos(Psi), 0],[sym.cos(Theta), 0, 1]])
B_inv = sym.simplify(B_a.inv())

Jq =  np.vstack((Jv_3,Jw_3))

X = sym.eye(3).col_join(sym.zeros(3,3))
X = X.row_join(sym.zeros(3,3).col_join(B_inv))
Ja = sym.simplify(X * Jq)

# inverse analytische Jacobimatrix --> Pseudoinverse
Ja = Ja.T

#Ja_t = sym.simplify(Ja.T*(Ja*Ja.T+lamda**2*sym.eye(3)).inv())
# Hilfskonstrukt zum Invertieren

Jat = Ja*Ja.T+lamda**2*sym.eye(3)
Jat_inv = Ainv.subs([(a11,Jat[0,0]),(a12,Jat[0,1]),(a13,Jat[0,2]), (a21,Jat[1,0]), (a22,Jat[1,1]), (a23,Jat[1,2]), (a31,Jat[2,0]), (a32,Jat[2,1]), (a33,Jat[2,2]) ])
Ja_t = Ja.T*Jat_inv


# zeitliche Ableitung der analytischen Jacobimatrix
Ja_d_q1 = sym.diff(Ja,q1)
Ja_d_q2 = sym.diff(Ja,q2)
Ja_d_q3 = sym.diff(Ja,q3)
Ja_diff = sym.simplify(Ja_d_q1*qd1 + Ja_d_q2*qd2 + Ja_d_q3)





######## ***************************************  
## 10.  Multivariable Control
######## ***************************************  

aq1,aq2,aq3= sym.symbols("aq1 aq2 aq3")
aq = sym.Matrix([aq1,aq2,aq3])

# Inverse Regelungsfunktion
#u_regler_ext = sym.simplify(M*aq+(C+B+R)*qd +gv.T)
u_regler_ext = M*aq+(C+B+R)*qd +gv.T