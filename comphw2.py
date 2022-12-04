import math as math

import matplotlib.pyplot as plt
import numpy as np

#For storing estimate values
TrapezoidalF1Vals = []
TrapezoidalF2Vals = []
TrapezoidalF3Vals = []

SimpsonsF1Vals = []
SimpsonsF2Vals = []
SimpsonsF3Vals = []

GaussF1Vals = []
GaussF2Vals = []
GaussF3Vals = []

#Absolute error for each method and function
TrapezoidalF1Errors = []
TrapezoidalF2Errors = []
TrapezoidalF3Errors = []

SimpsonsF1Errors = []
SimpsonsF2Errors = []
SimpsonsF3Errors = []

GaussF1Errors = []
GaussF2Errors = []
GaussF3Errors = []

#Rate of convergence for each method and function
ROCF1Trapezoidal = []
ROCF2Trapezoidal = []
ROCF3Trapezoidal = []

ROCF1Simpsons = []
ROCF2Simpsons = []
ROCF3Simpsons = []

ROCF1Gauss = []
ROCF2Gauss = []
ROCF3Gauss = []

#True integral values for each function obatined from WolfRamAlpha
TrueF1Val = -0.098363164308346596735 
TrueF2Val = 0.11033333333333333333
TrueF3Val = 1.7688741449316336902

#Def of each function
def f1(x):
    return np.cos(np.pi * x)

def f2(x):
    return np.piecewise(x, [x < 0, x >= 0], [lambda x: -(x * x), lambda x: x * x])

def f3(x):
    return np.exp((-(x**2)) / 2 )

#Trapezoidal def
def trapezoidal(f,n):
    h = (1.1 - -1) / n
    
    val = f(-1) + f(1.1)
    
    for i in range(1,n):
        k = -1 + i*h
        val += 2 * f(k)

    val *= h/2
    
    return val

#Simpsons def
def simpsons(f,n):
     
    h = (1.1 - -1) / n

    val = f(-1) + f(1.1)

    for i in range(1,n):
        k = -1 + i*h
        
        if i%2 == 0:
            val = val + (2 * f(k))

        elif i%2 != 0:
            val = val + (4 * f(k))
    
    val *= h/3

    
    return val

#Gauss n = 3 def

#Unable to solve for correct code
#def gauss(f,n):
#    A = []
#    x = []
#    a = -1
#    b = 1.1

#    #Xi and Ai values obtained from textbook pg 495
#    x0 = 0
#    x1 = 0.538469310105683 
#    x2 = 0.906179845938664
#    A0 = 0.568888888888889
#    A1 = 0.478628670499366
#    A2 = 0.236926885056189
#
#    x.append(x1)
#    x.append(x2)
#    A.append(A1)
#    A.append(A2)
#
#   h = (b - a) / n
#    S = 0
#
#    a = -1
#    b = 0.5
#    c = 0.5
#    d = 1.1
#
#   u = (a + b) /2
#   S = (A0 * f(u))
#
#    for i in range (0,2):
#            u = ((b-a) * x[i] + a + b) /2
#            v = ((a-b) * x[i] + a + b) / 2
##            S = S + (A[i]* (f(u) + f(v)))
 #           
 #   u = (c + d) /2
 #   Z = (A0 * f(u))
#
  #  for i in range (0,2):
  #          u = ((d-c) * x[i] + c + d) /2
  #          v = ((c-d) * x[i] + c + d) / 2
  #          Z = Z + (A[i]* (f(u) + f(v)))
#
#
   # print(S)
  #  print(Z)
#
  #  return S+Z

#Incorrect, but working Gaussian quad code

def gauss(f,n):
    half = float(1.1 - -1)/2.
    mid = (-1 + 1.1)/2.
    [x,w] = np.polynomial.legendre.leggauss(n)
    result =  0.0
    for i in range(n):
        result += w[i] * f(half*x[i] + mid)
    result *= half
    return result


for i in range(1,10):
    n = 2 ** i

    TrapezoidalF1Vals.append(trapezoidal(f1,n))
    TrapezoidalF2Vals.append(trapezoidal(f2,n))
    TrapezoidalF3Vals.append(trapezoidal(f3,n))

    SimpsonsF1Vals.append(simpsons(f1,n))
    SimpsonsF2Vals.append(simpsons(f2,n))
    SimpsonsF3Vals.append(simpsons(f3,n))

    GaussF1Vals.append(gauss(f1,n))
    GaussF2Vals.append(gauss(f2,n))
    GaussF3Vals.append(gauss(f3,n))


#manually checking the vals, for function 1, Gauss was closest to true value, Simpsons for function 2, and Gauss for function 3

#Set each corresponding end value to true value obtained from wolfram Alpha to avoid divide by 0 errors

GaussF1Vals[8] = -0.098363164308346596735
SimpsonsF2Vals[8] = 0.11033333333333333333
GaussF3Vals[8] = 1.7688741449316336902

#Getting errors for 2^i sub intervals up to i = 9; 512 subintervals
for i in range (1, 10):
    n = 2 ** i
    
    TrapezoidalF1Errors.append(abs(trapezoidal(f1,n) - GaussF1Vals[8]))
    TrapezoidalF2Errors.append(abs(trapezoidal(f2,n) - SimpsonsF2Vals[8]))
    TrapezoidalF3Errors.append(abs(trapezoidal(f3,n) - GaussF3Vals[8]))

    SimpsonsF1Errors.append(abs(simpsons(f1,n) - GaussF1Vals[8]))
    SimpsonsF2Errors.append(abs(simpsons(f2,n) - SimpsonsF2Vals[8]))
    SimpsonsF3Errors.append(abs(simpsons(f3,n) - GaussF3Vals[8]))

    GaussF1Errors.append(abs(gauss(f1,n) - GaussF1Vals[8]))
    GaussF2Errors.append(abs(gauss(f2,n) - SimpsonsF2Vals[8]))
    GaussF3Errors.append(abs(gauss(f3,n) - GaussF3Vals[8]))


#Calculating ROC for eachh method and function
for i in range (0,8):
    
    ROCF1Trapezoidal.append(math.log(TrapezoidalF1Errors[i] / TrapezoidalF1Errors[i+1]) / math.log(2))
    ROCF2Trapezoidal.append(math.log(TrapezoidalF2Errors[i] / TrapezoidalF2Errors[i+1]) / math.log(2))
    ROCF3Trapezoidal.append(math.log(TrapezoidalF3Errors[i] / TrapezoidalF3Errors[i+1]) / math.log(2))

    ROCF1Simpsons.append(math.log(SimpsonsF1Errors[i] / SimpsonsF1Errors[i+1]) / math.log(2))
    ROCF2Simpsons.append(math.log(SimpsonsF2Errors[i] / SimpsonsF2Errors[i+1]) / math.log(2))
    ROCF3Simpsons.append(math.log(SimpsonsF3Errors[i] / SimpsonsF3Errors[i+1]) / math.log(2))

    ROCF1Gauss.append(math.log(GaussF1Errors[i] / GaussF1Errors[i+1]) / math.log(2))
    ROCF2Gauss.append(math.log(GaussF2Errors[i] / GaussF2Errors[i+1]) / math.log(2))
    ROCF3Gauss.append(math.log(GaussF3Errors[i] / GaussF3Errors[i+1]) / math.log(2))


#Making CSV table

import csv
with open('MethodIntegrand.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Method', 'Trapezoidal', 'Trapezoidal', 'Trapezoidal', ' Simpsons', ' Simpsons', ' Simpsons', 'Gauss3', 'Gauss3', 'Gauss3'])
    writer.writerow(['Integrand', 'Consine', 'C1 function','normal dist', 'Consine', 'C1 function','normal dist', 'Consine', 'C1 function','normal dist'])
    writer.writerow([2, ROCF1Trapezoidal[0],ROCF2Trapezoidal[0],ROCF3Trapezoidal[0], ROCF1Simpsons[0],ROCF2Simpsons[0],ROCF3Simpsons[0], ROCF1Gauss[0], ROCF2Gauss[0], ROCF3Gauss[0]])
    writer.writerow([4,ROCF1Trapezoidal[1],ROCF2Trapezoidal[1],ROCF3Trapezoidal[1], ROCF1Simpsons[1],ROCF2Simpsons[1],ROCF3Simpsons[1], ROCF1Gauss[1], ROCF2Gauss[1], ROCF3Gauss[1]])
    writer.writerow([8,ROCF1Trapezoidal[2],ROCF2Trapezoidal[2],ROCF3Trapezoidal[2], ROCF1Simpsons[2],ROCF2Simpsons[2],ROCF3Simpsons[2], ROCF1Gauss[2], ROCF2Gauss[2], ROCF3Gauss[2]])
    writer.writerow([16,ROCF1Trapezoidal[3],ROCF2Trapezoidal[3],ROCF3Trapezoidal[3], ROCF1Simpsons[3],ROCF2Simpsons[3],ROCF3Simpsons[3], ROCF1Gauss[3], ROCF2Gauss[3], ROCF3Gauss[3]])
    writer.writerow([32,ROCF1Trapezoidal[4],ROCF2Trapezoidal[4],ROCF3Trapezoidal[4], ROCF1Simpsons[4],ROCF2Simpsons[4],ROCF3Simpsons[4], ROCF1Gauss[4], ROCF2Gauss[4], ROCF3Gauss[4]])
    writer.writerow([64,ROCF1Trapezoidal[5],ROCF2Trapezoidal[5],ROCF3Trapezoidal[5], ROCF1Simpsons[5],ROCF2Simpsons[5],ROCF3Simpsons[5], ROCF1Gauss[5], ROCF2Gauss[5], ROCF3Gauss[5]])
    writer.writerow([128,ROCF1Trapezoidal[6],ROCF2Trapezoidal[6],ROCF3Trapezoidal[6], ROCF1Simpsons[6],ROCF2Simpsons[6],ROCF3Simpsons[6], ROCF1Gauss[6], ROCF2Gauss[6], ROCF3Gauss[6]])
    writer.writerow([256,ROCF1Trapezoidal[7],ROCF2Trapezoidal[7],ROCF3Trapezoidal[7], ROCF1Simpsons[7],ROCF2Simpsons[7],ROCF3Simpsons[7], ROCF1Gauss[7], ROCF2Gauss[7], ROCF3Gauss[7]])


with open('IntegrandMethod.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Integrand', 'consine', 'consine', 'consine', 'C1 function', 'C1 function', 'C1 function', 'normal dist', 'normal dist', 'normal dist'])
    writer.writerow(['Method','Trapezoidal', 'Simpson', 'Gauss3','Trapezoidal', 'Simpson', 'Gauss3','Trapezoidal', 'Simpson', 'Gauss3'])
    writer.writerow([2, ROCF1Trapezoidal[0],ROCF1Simpsons[0], ROCF1Gauss[0],  ROCF2Trapezoidal[0],ROCF2Simpsons[0], ROCF2Gauss[0],ROCF3Trapezoidal[0],ROCF3Simpsons[0], ROCF3Gauss[0]])
    writer.writerow([4, ROCF1Trapezoidal[1],ROCF1Simpsons[1], ROCF1Gauss[1],  ROCF2Trapezoidal[1],ROCF2Simpsons[1], ROCF2Gauss[1],ROCF3Trapezoidal[1],ROCF3Simpsons[1], ROCF3Gauss[1]])
    writer.writerow([8, ROCF1Trapezoidal[2],ROCF1Simpsons[2], ROCF1Gauss[2],  ROCF2Trapezoidal[2],ROCF2Simpsons[2], ROCF2Gauss[2],ROCF3Trapezoidal[2],ROCF3Simpsons[2], ROCF3Gauss[2]])
    writer.writerow([16, ROCF1Trapezoidal[3],ROCF1Simpsons[3], ROCF1Gauss[3],  ROCF2Trapezoidal[3],ROCF2Simpsons[3], ROCF2Gauss[3],ROCF3Trapezoidal[3],ROCF3Simpsons[3], ROCF3Gauss[3]])
    writer.writerow([32, ROCF1Trapezoidal[4],ROCF1Simpsons[4], ROCF1Gauss[4],  ROCF2Trapezoidal[4],ROCF2Simpsons[4], ROCF2Gauss[4],ROCF3Trapezoidal[4],ROCF3Simpsons[4], ROCF3Gauss[4]])
    writer.writerow([64, ROCF1Trapezoidal[5],ROCF1Simpsons[5], ROCF1Gauss[5],  ROCF2Trapezoidal[5],ROCF2Simpsons[5], ROCF2Gauss[5],ROCF3Trapezoidal[5],ROCF3Simpsons[5], ROCF3Gauss[5]])
    writer.writerow([128, ROCF1Trapezoidal[6],ROCF1Simpsons[6], ROCF1Gauss[6],  ROCF2Trapezoidal[6],ROCF2Simpsons[6], ROCF2Gauss[6],ROCF3Trapezoidal[6],ROCF3Simpsons[6], ROCF3Gauss[6]])
    writer.writerow([256, ROCF1Trapezoidal[7],ROCF1Simpsons[7], ROCF1Gauss[7],  ROCF2Trapezoidal[7],ROCF2Simpsons[7], ROCF2Gauss[7],ROCF3Trapezoidal[7],ROCF3Simpsons[7], ROCF3Gauss[7]])
