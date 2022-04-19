# Sofia Palacios Cuevas
# Github: @spcJBS 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class LorenzAttractor:
    def Lorenz(self,t, X, sigma, beta, rho):
        'Lorenz eq'
        x,y,z = X
        dotX = sigma * (y - x)
        dotY = x * (rho - z) - y
        dotZ = x*y - beta*z
        return dotX, dotY, dotZ
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
        cmap = plt.cm.winter
        for i in range(0,n-s,s):
            ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=cmap(i/n), alpha=0.4)

        # Remove all the axis clutter, leaving just the curve.
        ax.set_axis_off()
        plt.savefig('myLorenz.png', dpi=DPI)
        plt.show()        
        
    def Rk4Lorenz():
        pass
    def EulerLorenz():
        pass



def main():
    tmax, n = 100, 10000
    x0, y0, z0 = 0,1,1.05
    sigma, beta, rho = 10, 8/3, 28
    WIDTH, HEIGHT, DPI = 1000, 750, 100
      
    la = LorenzAttractor()

    la.OdeLorenz(x0,y0,z0, tmax,n, sigma, beta, rho,WIDTH, HEIGHT, DPI)
    

if __name__ == "__main__":
    main()