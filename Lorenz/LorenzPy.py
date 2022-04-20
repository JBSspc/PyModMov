# Sofia Palacios Cuevas
# Github: @JBSspc 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from sympy.abc import x,y,z


class LorenzAttractor:

    # Lorenz attractor equations
    def Lorenz(self,t, X, sigma, beta, rho): 
        x,y,z = X
        dotX = sigma * (y - x)      # x' = sigma*(y-x)
        dotY = x * (rho - z) - y    # y' = x*(rho-z)-y
        dotZ = x*y - beta*z         # z' = x*y-beta*z
        return dotX, dotY, dotZ

    # Using scipy.integrate (solve_ivp) ==================================================================================================

    def OdeLorenz(self,x0,y0,z0, tmax,n, sigma, beta, rho,WIDTH, HEIGHT, DPI): 
        nSol = solve_ivp(self.Lorenz, (0,tmax), (x0,y0,z0), args = (sigma, beta, rho), dense_output=True)
        t = np.linspace(0, tmax,n)
        x,y,z = nSol.sol(t)
        # Plot the Lorenz attractor using a Matplotlib 3D projection.
        fig = plt.figure(facecolor='k', figsize=(WIDTH/DPI, HEIGHT/DPI))
        ax = fig.gca(projection='3d')
        ax.set_facecolor('k')
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

        # Make the line multi-coloured by plotting it in segments of length s which
        # change in colour across the whole time series.
        s = 10
        cmap = plt.cm.cool
        for i in range(0,n-s,s):
            ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=cmap(i/n), alpha=0.4)

        # Remove all the axis clutter, leaving just the curve.
        ax.set_axis_off()
        plt.savefig('myLorenz.png', dpi=DPI)
        plt.show()      
        
    # Runge-Kutta 4th order =============================================================================================================

    def L4rk(self,x,y,z, sigma, beta, rho):
        dotX = sigma * (y - x)      # x' = sigma*(y-x)
        dotY = x * (rho - z) - y    # y' = x*(rho-z)-y
        dotZ = x*y - beta*z         # z' = x*y-beta*z
        return dotX, dotY, dotZ
  
    def Euler(self,h,N,sigma, beta, rho):
        xe = np.empty((N+1,))
        ye = np.empty((N+1,))
        ze = np.empty((N+1,))

        xe[0], ye[0], ze[0] = (0.1,0.1,0.1)

        for i in range(N):
          x_dot, y_dot, z_dot = self.L4rk(xe[i], ye[i], ze[i],sigma, beta, rho)
          xe[i+1] = xe[i] + (x_dot*h)
          ye[i+1] = ye[i] + (y_dot*h)
          ze[i+1] = ze[i] + (z_dot*h)

        fig = plt.figure()
        ax = fig.gca(projection = '3d')

        ax.plot(xe, ye, ze)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Lorenz Attractor (Euler)')
        plt.savefig('LorenzEuler.png', dpi=100)
        plt.show()
        # 'https://www.youtube.com/watch?v=kAvJRF9BeiA'

    # Euler Method ======================================================================================================================
            
    def RK4(self, N,h, sigma, beta,rho):
        xrk4 = np.empty((N+1,)) # prealojamos
        yrk4 = np.empty((N+1,))
        zrk4 = np.empty((N+1,))

        xrk4[0], yrk4[0], zrk4[0] = (0.1,0.1,0.1)

        for i in range(N):
          x_dot, y_dot, z_dot = self.L4rk(xrk4[i],  yrk4[i], zrk4[i],sigma, beta, rho)
         
          k11 = x_dot
          k21 = x_dot + 0.5*h*k11
          k31 = x_dot + 0.5*h*k21
          k41 = x_dot +h*k31
          xrk4[i+1] = xrk4[i] + (h/6)*(k11 + 2*k21 + 2*k31 + k41)
          
          k12 = y_dot
          k22 = y_dot + 0.5*h*k12
          k32 = y_dot + 0.5*h*k22
          k42 = y_dot +h*k32
          yrk4[i+1] = yrk4[i] + (h/6)*(k12 + 2*k22 + 2*k32 + k42)
          
          k13 = z_dot
          k23 = z_dot + 0.5*h*k13
          k33 = z_dot + 0.5*h*k23
          k43 = z_dot +h*k33
          zrk4[i+1] = zrk4[i] + (h/6)*(k13 + 2*k23 + 2*k33 + k43)

        fig = plt.figure()
        ax = fig.gca(projection = '3d')

        ax.plot(xrk4,yrk4,zrk4)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Lorenz Attractor (RK4)')
        plt.savefig('LorenzRK4.png', dpi=100)
        plt.show()
        



def main():
    tfin = 200
    h = 0.01
    N = round((tfin-h)/h)
    tmax, n = 200, 10000
    x0, y0, z0 = 0.1,0.1,0.1
    sigma, beta, rho = 10, 8/3, 28
    WIDTH, HEIGHT, DPI = 1000, 750, 100
      
    la = LorenzAttractor()

    la.OdeLorenz(x0,y0,z0, tmax,n, sigma, beta, rho,WIDTH, HEIGHT, DPI)
    la.Euler(h,N,sigma, beta, rho)
    la.RK4(N,h, sigma, beta,rho)
    

if __name__ == "__main__":
    main()