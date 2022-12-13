## Titel:    Grundlagen Robot Control WS 2023
## Author:   Stefan Kaufmann
## Date:     17.11.2022

## Bib
import numpy as np
import sympy as  sym


######## ***************************************  
## 1. Transformation 2R Kette
######## ***************************************  
q1,q2,ai,a1,a2,di,d1,d2,alpha,alpha1,alpha2 = sym.symbols("q1 q2 ai a1 a2 di d1 d2 alpha alpha1 alpha2")

## For nice looking equations
sym.init_printing()


D = sym.Matrix([[1, 0, 0, 0],[ai, 1, 0, 0],[0, 0, sym.cos(alpha), -sym.sin(alpha)],[di, 0, sym.sin(alpha), sym.cos(alpha)]])


M1 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q1),-sym.sin(q1),0],[0,sym.sin(q1), sym.cos(q1), 0],[0, 0, 0, 1]])
G1 = sym.Matrix([[1, 0, 0, 0],[a1, 1, 0, 0],[0, 0, sym.cos(alpha1),-sym.sin(alpha1)],[d1, 0, sym.sin(alpha1), sym.cos(alpha1)]])
M2 = sym.Matrix([[1, 0, 0, 0],[0,  sym.cos(q2),-sym.sin(q2),0],[0,sym.sin(q2), sym.cos(q2), 0],[0, 0, 0, 1]])
G2 = sym.Matrix([[1, 0, 0, 0],[a2, 1, 0, 0],[0, 0, sym.cos(alpha2),-sym.sin(alpha2)],[d2, 0, sym.sin(alpha2), sym.cos(alpha2)]])
T02 = M1*G1*M2*G2
T01 = M1*G1


######## ***************************************  
## 2. Aufstellen der Jacobimatrix
######## ***************************************  


T2 = T02.subs({a1:3,a2:0,alpha1:sym.pi/6, alpha2:0,d1:1,d2:0,q1:sym.pi/3,q2:sym.pi/4})
x0 = sym.Matrix([1, 0, 0, 0])    # 1 ... ist für die Homogenisierung

x1 = T2*x0





