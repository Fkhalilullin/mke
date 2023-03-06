from math import *
import numpy as np
from numpy.linalg import inv
import sympy as sp
import pandas as pd

a = sp.symbols('a')
y = sp.symbols('y')
z = sp.symbols('z')

S=a**2

d = {
        'узел 1': [1,1,2,2,3,3,4,4,7,7,8,8,9,9,10,10], 
        'узел 2': [10,9,9,8,8,7,7,6,15,14,14,13,13,12,12,11],
        'узел 3': [9,2,8,3,7,4,6,5,6,15,7,14,8,13,9,12], 
        'y1': [0,0,0,0,0,0,0,0,a,a,a,a,a,a,a,a], 
        'y2': [a,a,a,a,a,a,a,a,2*a,2*a,2*a,2*a,2*a,2*a,2*a,2*a], 
        'y3': [a,0,a,0,a,0,a,0,a,2*a,a,2*a,a,2*a,a,2*a], 
        'z1': [0,0,a,a,2*a,2*a,3*a,3*a,3*a,3*a,2*a,2*a,a,a,0,0],
        'z2': [0,a,a,2*a,2*a,3*a,3*a,4*a,4*a,3*a,3*a,2*a,2*a,a,a,0],
        'z3': [a,a,2*a,2*a,3*a,3*a,4*a,4*a,4*a,4*a,3*a,3*a,2*a,2*a,a,a]
    }

elements = pd.DataFrame(data=d,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

elements = elements.assign(a1 = elements.y2*elements.z3-elements.y3*elements.z2)
elements = elements.assign(a2 = elements.y3*elements.z1-elements.y1*elements.z3)
elements = elements.assign(a3 = elements.y1*elements.z2-elements.y2*elements.z1)

elements = elements.assign(b1 = elements.z2-elements.z3)
elements = elements.assign(b2 = elements.z3-elements.z1)
elements = elements.assign(b3 = elements.z1-elements.z2)

elements = elements.assign(c1 = elements.y3-elements.y2)
elements = elements.assign(c2 = elements.y1-elements.y3)
elements = elements.assign(c3 = elements.y2-elements.y1)

elements = elements.assign(L1 = 1/(2*S)*(elements.a1+elements.b1*y+elements.c1*z))
elements = elements.assign(L2 = 1/(2*S)*(elements.a2+elements.b2*y+elements.c2*z))
elements = elements.assign(L3 = 1/(2*S)*(elements.a3+elements.b3*y+elements.c3*z))

K=[]

for i in range(1,len(elements)+1):
    e1=np.matrix([
        [elements.b1[i]*elements.b1[i],elements.b1[i]*elements.b2[i],elements.b1[i]*elements.b3[i]],
        [elements.b2[i]*elements.b1[i],elements.b2[i]*elements.b2[i],elements.b2[i]*elements.b3[i]],
        [elements.b3[i]*elements.b1[i],elements.b3[i]*elements.b2[i],elements.b3[i]*elements.b3[i]]])
    e2=np.matrix([
        [elements.c1[i]*elements.c1[i],elements.c1[i]*elements.c2[i],elements.c1[i]*elements.c3[i]],
        [elements.c2[i]*elements.c1[i],elements.c2[i]*elements.c2[i],elements.c2[i]*elements.c3[i]],
        [elements.c3[i]*elements.c1[i],elements.c3[i]*elements.c2[i],elements.c3[i]*elements.c3[i]]])
    e=e1+e2
    K.append(1/(2*S)*e)
elements = elements.assign(Ki=K)

elements['P'] = 1/3*a**2*np.ones(16)

f =[1,2,3,4,7,8,9,10]

K_ob = np.zeros((8,8))

p1 = np.zeros(8)
#print(K_ob)
for str in range(1,len(elements)+1):
    str_i =np.array(np.array(elements.iloc[str-1:str,0:3])[0])
    k=0
    a=-1
    b=-1
    c=-1
    print(str_i)
#     print(f)
    for i in f:
        for j in range(0,3):
            #print(str_i)
            #print(str_i[j])
            if str_i[j] == i:
                K_ob[f.index(i),f.index(i)]+=K[str-1][j,j]
#               print(f.index(i))
                if j==0:a=i
                if j==1:b=i
                if j==2:c=i
                k+=1
                p1[f.index(i)]+=1
    if k>1:
        str_i = list(str_i)
        if(a!=-1 and b!=-1):
            K_ob[f.index(a),f.index(b)]+=K[str-1][str_i.index(a),str_i.index(b)]
            K_ob[f.index(b),f.index(a)]+=K[str-1][str_i.index(b),str_i.index(a)]
        if(a!=-1 and c!=-1):
            K_ob[f.index(a),f.index(c)]+=K[str-1][str_i.index(a),str_i.index(c)]
            K_ob[f.index(c),f.index(a)]+=K[str-1][str_i.index(c),str_i.index(a)]
        if(b!=-1 and c!=-1):
            K_ob[f.index(b),f.index(c)]+=K[str-1][str_i.index(b),str_i.index(c)]
            K_ob[f.index(c),f.index(b)]+=K[str-1][str_i.index(c),str_i.index(b)]
   
# print(k)
# print(a,b,c)
# print(K[str-1])
# print(K_ob)

print(p1)
a = sp.symbols('a')
P = 1/3*p1

F = np.linalg.solve(K_ob, P)
J=0
for i in range(0,8):
    J+=1/2*4*2*1/3*a**4*F[i]*p1[i]

print(J)