# Sofia Palacios Cuevas
# Github: @spcJBS 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class LorenzAttractor:
    def OdeLorenz(self,x0,y0,z0, tmax,n, sigma, beta, rho):
        nSol = solve_ivp(Lorenz, (0,tmax), (x0,y0,z0), args = (sigma, beta, rho), dense_output=True)
        t = np.linspace(0, tmax,n)
        return nSol
        
    def Rk4Lorenz():
        pass
    def EulerLorenz():
        pass

def Lorenz(t, X, sigma, beta, rho):
    'Lorenz eq'
    x,y,z = X
    dotX = sigma * (y - x)
    dotY = x * (rho - z) - y
    dotZ = x*y - beta*z
    return dotX, dotY, dotZ

def main():
    tmax, n = 100, 10000
    x0, y0, z0 = 0,1,1.05
    sigma, beta, rho = 10, 8/3, 28
    
    la = LorenzAttractor(x0,y0,z0, tmax,n, sigma, beta, rho)
    la.OdeLorenz()