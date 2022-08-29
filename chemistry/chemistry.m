% Daniel R. Reynolds
% SMU Mathematics
% Math 4370 / 6370
%-----------------------------------------------------------------
% Description:
%    Computes the equilibrium chemical densities at a number of
%    spatial locations, given a (random) background temperature
%    field.  The chemical rate equations and solution strategy
%    are in the subroutine chem_solver, which is called at every
%    spatial location.
%=================================================================

% clear all old variables from memory
clear

% declare variables to be used
n = 0;
i = 0;
maxit = 0;
its = 0;
lam = 0;
err = 0;
res = 0;
stime = 0;
ftime = 0;

% set solver input parameters
maxit = 1000000;
lam = 1.d-2;
err = 1.d-10;

% input the number of intervals
n = input('Enter the number of intervals (0 quits): ');
if (n < 1)
  return
end

% allocate temperature and solution arrays
T = zeros(n,1);
u = zeros(n,1);
v = zeros(n,1);
w = zeros(n,1);

% set random temperature field, initial guesses at chemical densities
T = rand(size(T));
u = 0.35*ones(size(u));
v = 0.1*ones(size(v));
w = 0.5*ones(size(w));

% start timer
stime = cputime;

% call solver over n intervals
for i=1:n
  [u(i),v(i),w(i),its,res] = chem_solver(T(i),u(i),v(i),w(i),lam,err,maxit);
  if (res < err)
    disp(sprintf('  i = %i,  its = %i',i,its));
  else
    disp(sprintf('  error: i = %i, its = %i, res = %g, u = %g, v = %g, w = %g',...
	               res,u(i),v(i),w(i)));
    return
  end
end

% stop timer
ftime = cputime;

% output solution time
disp(sprintf('     runtime = %g',ftime-stime));


% end program
