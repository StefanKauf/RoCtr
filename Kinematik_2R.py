## Titel:    Grundlagen Robot Control WS 2023
## Author:   Stefan Kaufmann
## Date:     17.11.2022

## Bib
import numpy as np
import sympy as  sym
import Parameter as param



######## ***************************************  
## 1. Transformation 2R Kette
######## ***************************************  

q1,q2,ai,a1,a2,di,d1,d2,alpha,alpha1,alpha2,l_s1,l_s2,l1,l2= sym.symbols("q1 q2 ai a1 a2 di d1 d2 alpha alpha1 alpha2 l_s1 l_s2 l1 l2")

## For nice looking equations
#sym.init_printing()

 
D = sym.Matrix([[1, 0, 0, 0],[ai, 1, 0, 0],[0, 0, sym.cos(alpha), -sym.sin(alpha)],[di, 0, sym.sin(alpha), sym.cos(alpha)]])
q = sym.Matrix([q1,q2])

M1 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q1),-sym.sin(q1),0],[0,sym.sin(q1), sym.cos(q1), 0],[0, 0, 0, 1]])
G1 = sym.Matrix([[1, 0, 0, 0],[a1, 1, 0, 0],[0, 0, sym.cos(alpha1),-sym.sin(alpha1)],[d1, 0, sym.sin(alpha1), sym.cos(alpha1)]])
M2 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q2),-sym.sin(q2),0],[0,sym.sin(q2), sym.cos(q2), 0],[0, 0, 0, 1]])
G2 = sym.Matrix([[1, 0, 0, 0],[a2, 1, 0, 0],[0, 0, sym.cos(alpha2),-sym.sin(alpha2)],[d2, 0, sym.sin(alpha2), sym.cos(alpha2)]])
T02 = M1*G1*M2*G2
T01 = M1*G1




######## ***************************************  
## 2. Aufstellen der Jacobimatrix
######## ***************************************  

# Geschwindigkeit in x,y Richtung
Jv_01 = T01[1:3,0].jacobian(q)
Jv_02 = T02[1:3,0].jacobian(q)    # 1:3  --> Jacobimatrix bis zur 2ten Dimension, 1 Zeile ist die homogenisierung


# Rotationsgeshwindigkeit

Jw_01 = M1*sym.Matrix([0, 0, 0, 1])+M1[:,0]
Jw_02 = Jw_01.row_join(T01*sym.Matrix([0, 0, 0, 1]))
Jw_01 = Jw_01.row_join(sym.Matrix([0,0,0,0]))

# erste Zeile Löschen
Jw_01 = Jw_01[1:4,:]
Jw_02 = Jw_02[1:4,:]





# Substitueren mit den Parametern

Jv_1 = Jv_01.subs({a1:l_s1,alpha1:0,d1:0})
Jv_2 = Jv_02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})

Jw_1 = Jw_01.subs({a1:l_s1,alpha1:0,d1:0})
Jw_2 = Jw_02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})




######## ***************************************  
## 3. Aufstellen der D-Matrix
######## ***************************************  

m1,m2,I1,I2,J1,J2,B1,B2,R1,R2,km1,km2,kb1,kb2,r1,r2 =sym.symbols("m1 m2 I1 I2 J1 J2 B1 B2 R1 R2 km1 km2 kb1 kb2 r1 r2")

D = m1*Jv_1.T*Jv_1 + Jw_1.T*I1*Jw_1 + m2*Jv_2.T*Jv_2 + Jw_2.T*I2*Jw_2 
#sym.simplify(D)

######## ***************************************  
## 4. Aufstellen der C-Matrix
######## ***************************************  

qd1, qd2, qdd1, qdd2 = sym.symbols("qd1 qd2 qdd1 qdd2")
qd = sym.Matrix([qd1,qd2])
qdd = sym.Matrix([qdd1,qdd2])


c = sym.MutableDenseNDimArray(np.zeros((2,)*3))

n = 2
for i in range(n):
    for j in range(n):
        for k in range(n):            
            c[i,j,k] = 1/2*(sym.diff(D[k,j],q[i])+ sym.diff(D[k,i],q[j]) - sym.diff(D[i,j],q[k]) )
           


C = sym.zeros(n,n)

for i in range(n):
    for k in range(n):
        for j in range(n):
            C[k,i] = C[k,i] + c[i,j,k]*qd[j]
        
    
#sym.simplify(C)

######## ***************************************  
## 5. Aufstellen der g-Vektors  --> potentielle Energie
######## ***************************************  
g = sym.symbols("g")

T01s = T01.subs({a1:l_s1,alpha1:0,d1:0})
T02s = T02.subs({a1:l1,a2:l_s2,alpha1:0, alpha2:0,d1:0,d2:0})

# Potentielle Energie
P = sym.Matrix([0, 0, g, 0]).T*T01s[:,0]*m1 + sym.Matrix([0, 0, g, 0]).T*T02s[:,0]*m2

gv = P.jacobian(q)



######## ***************************************  
## 6. Bewegungsgleichung einfaches Modell
######## ***************************************  
u1,u2 = sym.symbols("u1 u2")
tau = sym.Matrix([u1,u2])
#tau = D*qdd + C*qd + gv
'''
qdd_Modell = sym.Function('qdd_modell')(qd1,qd2,q1,q2,u1,u2)
qdd_modell = sym.simplify(D.inv()*(tau-C*qd-gv.T))
'''

######## ***************************************  
## 7. Erweitertes Modell
######## ***************************************  

J = sym.Matrix([[J1*r1**2, 0],[0, J2*r2**2]])
B = sym.Matrix([[r1**2*B1, 0],[0, B2*r2**2]])
R = sym.Matrix([[km1*kb1/R1, 0],[0, km2*kb2/R2]])

#tau = sym.Matrix([r1*km1/R1*u1, r2*km2/R2*u2])

# tau = (D+J)qdd + (C+B+R)qd +gv
M = D+J

qdd_ext = M.inv()*(tau-(C+B+R)*qd - gv.T)


######## ***************************************  
## 8.  Einsetzen der Parameter
######## ***************************************  

'''

qdd_modell_subs = sym.simplify(qdd_modell.subs({l1:param.l1,l2:param.l2,I1:param.I1,I2:param.I2,m1:param.m1,m2:param.m2,l_s1:param.l_s1,l_s2:param.l_s2,g:param.g,J1:param.J1,J2:param.J2, B1:param.B1, B2:param.B2,R1:param.R1,R2:param.R2,r1:param.r1,r2:param.r2,km1:param.km1,km2:param.km2,kb1:param.kb1,kb2:param.kb2}))
f_modell = sym.lambdify([qd1, qd2,q1,q2,u1,u2], qdd_modell_subs)

qdd_ext_subs = sym.simplify(qdd_ext.subs({l1:param.l1,l2:param.l2,I1:param.I1,I2:param.I2,m1:param.m1,m2:param.m2,l_s1:param.l_s1,l_s2:param.l_s2,g:param.g,J1:param.J1,J2:param.J2, B1:param.B1, B2:param.B2,R1:param.R1,R2:param.R2,r1:param.r1,r2:param.r2,km1:param.km1,km2:param.km2,kb1:param.kb1,kb2:param.kb2}))
f_modell_ext = sym.lambdify([qd1, qd2,q1,q2,u1,u2], qdd_ext_subs)

'''




######## ***************************************  
## 9.  analytische Jacobimatrix
######## ***************************************  
Psi,Theta = sym.symbols("Psi Theta")
B_a = sym.Matrix([[sym.cos(Psi)*sym.cos(Theta),-sym.sin(Theta), 0],[sym.sin(Psi)*sym.sin(Theta), sym.cos(Psi), 0],[sym.cos(Theta), 0, 1]])

J = Jv_2.col_join(sym.zeros(1,2)).col_join(Jw_2)

X = sym.eye(3).col_join(sym.zeros(3,3))
X = X.row_join(sym.zeros(3,3).col_join(B_a))
Ja = sym.simplify(X * J)

# inverse analytische Jacobimatrix --> Pseudoinverse
xd, yd, zd, xw, yw, zw, lamda= sym.symbols("x_dot y_dot z_dot x_omega y_omega z_omega lamda")
xe = sym.Matrix([xd, yd, zd, xw, yw, zw ])

Ja_t = Ja.T*(Ja*Ja.T+lamda**2*sym.eye(2))

# zeitliche Ableitung der analytischen Jacobimatrix

######## ***************************************  
## 9.  analytische Jacobimatrix
######## ***************************************  
'''
# Analytic Jacobi matrix
Ja = sym.simplify(J[0:2,:])


# Inverse Jacobian
Ja_inv = sym.simplify(Ja.inv())

# Time devireative of the jacobian
Ja_d_q1 = sym.diff(Ja,q1)
Ja_d_q2 = sym.diff(Ja,q2)
Ja_diff = sym.simplify(Ja_d_q1*qd1 + Ja_d_q2*qd2)

'''




######## ***************************************  
## 10.  Multivariable Control
######## ***************************************  

aq1,aq2 = sym.symbols("aq1 aq2")
aq = sym.Matrix([aq1,aq2])

# Inverse Regelungsfunktion
u_regler = sym.simplify(D*aq+C*qd +gv.T)
u_regler_ext = sym.simplify(M*aq+(C+B+R)*qd +gv.T)

#### Kartesisch in Gelenksraum
# analytische Jacobimatrix




