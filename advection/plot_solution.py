#!/usr/bin/env python3
#
# Inputs the file u_sol.txt to determine problem/parallelism information.
# Inputs the files u_sol.000.000, u_sol.001.000, ..., to read the output from
# each MPI process at each time step.  Concatenates these solutions together
# and displays plots for each solution.
#
# Implemented to run in parallel, spawning separate threads for each snapshot.
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 4370 / 6370

# imports
from multiprocessing import Process
from pylab import *
import numpy as np
from os import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

# set the desired output file type
filetype = '.png'


# set up some helper functions

###
def load_info():
    """Returns the mesh size, parallelism information, and total
       number of output times from the input file 'u_sol.txt':
          nx,ny,px,py,nt = load_info() """
    data = np.loadtxt("u_sol.txt", dtype=int)
    return data

###
def load_snapshot(tstep):
    """Returns the solution over the entire (parallel) mesh for
       a given time snapshot:
          t,u = load_snapshot(tstep) """

    # load the parallelism information
    nx,ny,px,py,nt = load_info()

    # check that tstep is allowed
    if (tstep > nt):
        print('load_snapshot error: tstep exceeds outputs!')
        return

    # allocate data for snapshot
    u = np.zeros((nx,ny), dtype=double)

    # determine local subdomain sizes
    nloc = nx//px   # // is like /, but enforces "flooring" of the result
    mloc = ny//py

    # iterate over processors, loading data into output
    for iproc in range(px*py):

        # determine data file name
        sproc = repr(iproc).zfill(3)
        stime = repr(tstep).zfill(3)
        outfile = 'u_sol.' + sproc + '.' + stime

        # load file, split into relevant parts
        data = np.loadtxt(outfile, dtype=double)
        nxloc = int(data[0])
        nyloc = int(data[1])
        pxloc = int(data[2])
        pyloc = int(data[3])
        t     = data[4]
        uloc = data[5:]
        uloc = np.reshape(uloc, (nxloc,nyloc), order='F')

        # insert data into output
        istart = pxloc*nloc
        iend   = istart+nxloc
        jstart = pyloc*mloc
        jend   = jstart+nyloc
        u[istart:iend, jstart:jend] = uloc

    return [t,u]

###
def plot_snapshot(tstep):
    """Plots the snapshot as both a 3D surface and as a filled contour plot """

    nx,ny,px,py,nt = load_info()

    # input solution at this time
    t,u = load_snapshot(tstep)

    # set string constants for output plots, current time, mesh size
    pname = 'sol.' + repr(tstep).zfill(3) + filetype
    cname = 'contour.' + repr(tstep).zfill(3) + filetype
    tstr = repr(round(t,4))
    nxstr = repr(nx)
    nystr = repr(ny)

    # set x and y meshgrid objects
    xspan = np.linspace(0.0, 1.0, nx)
    yspan = np.linspace(0.0, 1.0, ny)
    X,Y = np.meshgrid(xspan,yspan)

    # plot current solution as a surface, and save to disk
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.jet,
                    linewidth=0, antialiased=True, shade=True)
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    title('u(x,y) at t = ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    savefig(pname)
    plt.close()

    # plot current solution as a contour plot, and save to disk
    figure()
    h = imshow(u, origin='lower')
    colorbar(h, orientation='vertical')
    title('u(x,y) at t = ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    savefig(cname)
    plt.close()





##################

# here is the actual script

# input general problem information
nx,ny,px,py,nt = load_info()

# iterate over time steps
p = []
for it in range(nt+1):

    # launch a thread to create snapshot for this step
    #p = Process(target=plot_snapshot, args=(it,))
    #p.start()
    plot_snapshot(it)

# end of script
