#!/usr/bin/env python
#
# Testing script for running hybrid MPI+OpenMP tests on Xeon Phi
# 
# Daniel R. Reynolds
# SMU Mathematics
# Math 4370/6370
# 22 April 2017

# imports
import os
import subprocess
from time import time
import numpy as np
import pylab as plt


# testing parameters
MaxP = 256
Nint = 50000
MPIEXEC = '/opt/intel/compilers_and_libraries_2017.0.098/linux/mpi/intel64/bin/mpiexec'
EXEC = 'chemistry_hybrid_cpp.exe'
my_env = os.environ.copy()

# arrays for timing results
Procs = np.linspace(1,MaxP,MaxP)
MPITimes = np.zeros(MaxP)
OMPTimes = np.zeros(MaxP)
Hybrid2Times = np.zeros(MaxP//2)
Hybrid2Procs = np.zeros(MaxP//2)
Hybrid4Times = np.zeros(MaxP//4)
Hybrid4Procs = np.zeros(MaxP//4)
Hybrid8Times = np.zeros(MaxP//8)
Hybrid8Procs = np.zeros(MaxP//8)
Hybrid16Times = np.zeros(MaxP//16)
Hybrid16Procs = np.zeros(MaxP//16)


# run MPI-only tests
t = 1
for p in range(1,MaxP+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   MPITimes[p-1] = tend-tstart


# run OpenMP-only tests
p = 1
for t in range(1,MaxP+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   OMPTimes[t-1] = tend-tstart


# run Hybrid2 tests
p = 2
for t in range(1,MaxP//p+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   Hybrid2Times[t-1] = tend-tstart
   Hybrid2Procs[t-1] = t*p


# run Hybrid4 tests
p = 4
for t in range(1,MaxP//p+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   Hybrid4Times[t-1] = tend-tstart
   Hybrid4Procs[t-1] = t*p


# run Hybrid8 tests
p = 8
for t in range(1,MaxP//p+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   Hybrid8Times[t-1] = tend-tstart
   Hybrid8Procs[t-1] = t*p


# run Hybrid16 tests
p = 16
for t in range(1,MaxP//p+1):
   
   print '  '
   my_env["OMP_NUM_THREADS"] = str(t)
   tstart = time()
   subprocess.call([MPIEXEC, "-n", str(p), EXEC, str(Nint)], env=my_env)
   tend = time()
   Hybrid16Times[t-1] = tend-tstart
   Hybrid16Procs[t-1] = t*p


# define plotting function
def speedup_efficiency(procs,times,ftitle,fname):
   tserial = times[0]*procs[0]
   speedup = tserial/times
   efficiency = speedup/procs
   fig, ax1 = plt.subplots()
   ax1.plot(procs, speedup, 'b-')
   ax1.set_xlabel('procs')
   ax1.set_ylabel('speedup', color='b')
   ax1.tick_params('y', colors='b')
   ax2 = ax1.twinx()
   ax2.plot(procs, efficiency, 'r.')
   ax2.set_ylabel('efficiency', color='r')
   ax2.tick_params('y', colors='r')
   fig.tight_layout()
   plt.title(ftitle)
   plt.savefig(fname)
   plt.close()
   


# generate/save speedup/efficiency plots
speedup_efficiency(Procs,MPITimes,'MPI-Only Performance','mpi.png')
speedup_efficiency(Procs,OMPTimes,'OpenMP-Only Performance','omp.png')
speedup_efficiency(Hybrid2Procs,Hybrid2Times,'Hybrid-2 Performance','hybrid2.png')
speedup_efficiency(Hybrid4Procs,Hybrid4Times,'Hybrid-4 Performance','hybrid4.png')
speedup_efficiency(Hybrid8Procs,Hybrid8Times,'Hybrid-8 Performance','hybrid8.png')
speedup_efficiency(Hybrid16Procs,Hybrid16Times,'Hybrid-16 Performance','hybrid16.png')


# save data
np.savetxt('procs.txt', Procs);
np.savetxt('times_MPI.txt', MPITimes);
np.savetxt('times_OMP.txt', OMPTimes);
np.savetxt('procs_hybrid2.txt', Hybrid2Procs);
np.savetxt('times_hybrid2.txt', Hybrid2Times);
np.savetxt('procs_hybrid4.txt', Hybrid4Procs);
np.savetxt('times_hybrid4.txt', Hybrid4Times);
np.savetxt('procs_hybrid8.txt', Hybrid8Procs);
np.savetxt('times_hybrid8.txt', Hybrid8Times);
np.savetxt('procs_hybrid16.txt', Hybrid16Procs);
np.savetxt('times_hybrid16.txt', Hybrid16Times);


