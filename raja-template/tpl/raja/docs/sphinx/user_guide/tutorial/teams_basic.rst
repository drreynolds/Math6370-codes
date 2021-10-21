.. ##
.. ## Copyright (c) 2016-20, Lawrence Livermore National Security, LLC
.. ## and RAJA project contributors. See the RAJA/COPYRIGHT file
.. ## for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)
.. ##

.. _teamsbasic-label:

------------------------------
Team based loops (RAJA Teams)
------------------------------

Key RAJA features shown in the following examples:

  * ``RAJA::expt::launch`` method to create a run-time
    selectable host/device execution space.
  * ``RAJA::expt::loop`` methods to express algorithms
    in terms of nested for loops. 

In this example, we introduce the RAJA Teams framework and discuss
hierarchical loop-based parallelism. Development with RAJA Teams occurs
inside an execution space. The execution space is launched using the 
``RAJA::expt::launch`` method::

  RAJA::expt::launch<launch_policy>(RAJA::expt::ExecPlace ,
  RAJA::expt::Grid(RAJA::expt::Teams(Nteams,Nteams),
                        RAJA::expt::Threads(Nthreads,Nthreads)),
  [=] RAJA_HOST_DEVICE (RAJA::expt::LaunchContext ctx) {

    /* Express code here */

  });

The ``RAJA::expt::launch`` method is templated on both a host and a device launch policy.
As an example, the following constructs an execution space for a sequential 
and CUDA kernel::

  using launch_policy = RAJA::expt::LaunchPolicy
    <RAJA::expt::seq_launch_t, RAJA::expt::cuda_launch_t<false>>;

Kernel execution on either the host or device is driven by the first argument of
the method which takes a ``RAJA::expt::ExecPlace`` enum type, either ``HOST`` or ``DEVICE``. 
Similar to thread, and block programming models, RAJA Teams carries out
computation in a predefined compute grid made up of threads which are
then grouped into teams. The execution space is then enclosed by a host/device
lambda which takes a ``RAJA::expt::LaunchContext`` object. The ``RAJA::expt::LaunchContext``
may then be used to control the flow within the kernel, for example creating thread-team
synchronization points. 

Inside the execution space the ``RAJA::expt::loop`` methods enable developers
to express their code in terms of nested loops. The manner in which the loops
are executed depends on the template. Following the CUDA/HIP programming models
we follow a hierarchical structure in which outer loops are executed by thread-teams
and inner loops are executed by a thread in a team. 

.. literalinclude:: ../../../../examples/tut_teams_basic.cpp
   :start-after: // _team_loops_start
   :end-before: // _team_loops_end
   :language: C++
  
The mapping between the thread and teams to programming model depends on 
how they are defined. For example, we may define host and device mapping 
strategies as the following::

  using teams_x = RAJA::expt::LoopPolicy<RAJA::loop_exec,
                                         RAJA::cuda_block_x_direct>;
  using thread_x = RAJA::expt::LoopPolicy<RAJA::loop_exec, 
                                          RAJA::cuda_block_x_direct>;

In the example above the ``RAJA::expt::LoopPolicy`` struct holds both the host and
device loop mapping strategies. On the host, both the team/thread strategies expand
out to standard C-style loops for execution:

.. literalinclude:: ../../../../examples/tut_teams_basic.cpp
   :start-after: // _c_style_loops_start
   :end-before: // _c_style_loops_end
   :language: C++
   
On the device the ``teams_x/y`` policies will map loop iterations directly to 
CUDA thread blocks, while the ``thread_x/y`` policies will map loop iterations
directly to threads in a CUDA block. The CUDA equivalent is illustrated below:   

.. literalinclude:: ../../../../examples/tut_teams_basic.cpp
   :start-after: // _device_loop_start
   :end-before: // _device_loop_end
   :language: C++
   
The file RAJA/examples/tut_teams_basic.cpp contains the complete working example code.
