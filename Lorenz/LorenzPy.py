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

    # using scipy.integrate (solve_ivp) 
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
        
    # Runge-Kutta 4th order

    def L4rk(self,x,y,z, sigma, beta, rho):
        dotX = sigma * (y - x)      # x' = sigma*(y-x)
        dotY = x * (rho - z) - y    # y' = x*(rho-z)-y
        dotZ = x*y - beta*z         # z' = x*y-beta*z
        return dotX, dotY, dotZ
  
    def Rk4(self,h,N,sigma, beta, rho):
        xs = np.empty((N+1,))
        ys = np.empty((N+1,))
        zs = np.empty((N+1,))

        xs[0], ys[0], zs[0] = (0.1,0.1,0.1)

        for i in range(N):
          x_dot, y_dot, z_dot = self.L4rk(xs[i], ys[i], zs[i],sigma, beta, rho)
          xs[i+1] = xs[i] + (x_dot*h)
          ys[i+1] = ys[i] + (y_dot*h)
          zs[i+1] = zs[i] + (z_dot*h)

        fig = plt.figure()
        ax = fig.gca(projection = '3d')

        ax.plot(xs, ys, zs)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Lorenz Attractor (RK4)')
        plt.savefig('LorenzRK4.png', dpi=100)
        plt.show()

            
    def EulerLorenz():
        pass



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
    la.Rk4(h,N,sigma, beta, rho)
    

if __name__ == "__main__":
    main()