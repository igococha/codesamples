__author__ = 'Igor'

import numpy as np
import matplotlib.pyplot as plt

### One dimensional brownian motion (Wiener process)

def generate_brownian(T,n):
    # T =1; N=500
    # break T into N intervals
    # dt = spacing between samples
    t,dt = np.linspace(0,T,n+1,retstep=True)
    dW = np.random.normal(0.0,np.sqrt(dt),n+1)
    dW[0] = 0.0
    W = np.cumsum(dW)
    return t,W

def plot_brownian(t,W):
    plt.plot(t,W)
    plt.xlabel('t')
    plt.ylabel('W(t)')
    plt.title('Sample Wiener Process',weight='bold',size=16)

def simulate_brownian(T,n):
    t,W = generate_brownian(T,n)
    plot_brownian(t,W)

# phi(x) = x + (1/10)*sin(3PIx)
def plot_phi():
    n = 10
    nf = 1000

    # for the red points
    x = np.linspace(0,1,n)
    # for the curve
    xf = np.linspace(0,1,nf)

    phi = lambda x: x+0.1*np.sin(3*np.pi*x)

    y = phi(x)

    # default style 'b-' solid blue line
    plt.plot(xf,phi(xf),lw=2) # lw line width
    # 'ro' plot points as red circles
    plt.plot(x,y,'ro',ms=8) # ms marker size

    # The 'bases' (blobs) of the vlines and hlines
    plt.plot(x,[0]*n,'ko',ms=8)
    plt.plot([0]*n,y,'ko',ms=8)
    # Dotted vlines and hlines
    plt.vlines(x,[0],y,linestyle='--')
    plt.hlines(y,[0],x,linestyle='--')
    plt.title(r'$y=\phi(x)$',fontsize=24)
    plt.xlabel(r'$x$',fontsize=24)
    plt.ylabel(r'$y$',fontsize=24)
    plt.show()