function [u,v,w,its,res] = chem_solver(T,u,v,w,lam,err,maxit)
% Usage: [u,v,w,its,res] = chem_solver(T,u,v,w,lam,err,maxit)
%
%-----------------------------------------------------------------
% Daniel R. Reynolds
% SMU, Mathematics
% Math 4370/6370
% 7 February 2015
%-----------------------------------------------------------------
% Description: 
%    Computes the equilibrium chemical concentrations, given a 
%    background temperature field, of a simple reaction network:
%        u + x -> v  with rate k1
%        u + x <- v  with rate k2
%        v -> w      with rate k3
%        v <- w      with rate k4,
%    where we have assumed that the total concentration 
%    (x+u+v+w) = 1, and where k1(T), k2(T), k3(T), and k4(T) are 
%    the temperature-dependent coefficient functions,
%        k1(T) = exp(-5*T),
%        k2(T) = atan(5*(T-1/2))/3 + 1/2,
%        k3(T) = 1/cosh(5*(T-1/2)),
%        k4(T) = tanh(5*(T-1/2)^2).
%    Using the conservation relation, we write the constrained 
%    ODE system
%          x = 1 - u - v - w,
%        u_t = k2*v - k1*u*x,
%        v_t = k1*u*x - (k2+k3)*v + k4*w,
%        w_t = k3*v - k4*w.
%    Inserting the constraint equation into the rate equations, 
%    and setting the time derivatives equal to zero (to find 
%    equilibrium concentrations), we have the system
%        0 = k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
%        0 = k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
%        0 = k3*v - k4*w                        = fw(u,v,w),
%    where each of the rate coefficients are frozen at the fixed 
%    temperature T.  
%
%    To solve this system, we call a simple damped fixed-point 
%    iteration: given an initial guess X0, compute iterates
%        Xn = X{n-1} + lambda*f, 
%    where 0 < lambda <= 1 is the damping parameter.  We compute 
%    these iterates Xn until |f(Xn)| < err.
% 
% Arguments:
%    T - (input), temperature
%    u - (in/out), concentration (in: guess, out: solution)
%    v - (in/out), concentration (in: guess, out: solution)
%    w - (in/out), concentration (in: guess, out: solution)
%  lam - (in), damping parameter (lambda)
%  err - (in), nonlinear solver tolerance
%  maxit - (in), maximum allowed iterations
%  its - (out), # of iterations required for convergence
%  res - (out), final nonlinear residual (max norm)
%=================================================================


% declare local variables for this function
k1 = 0;
k2 = 0;
k3 = 0;
k4 = 0;
f = zeros(3,1);

% compute chemical rate coefficients
k1 = exp(-5*T);
k2 = atan(5*(T-0.5))/3 + 0.5;
k3 = 1/cosh(5*(T-0.5));
k4 = tanh(5*(T-0.5)^2);

% compute initial residual function
f = chem_residual;
res = max(abs(f));

% loop
for its = 0:maxit
   
   % check for convergence
%   disp(sprintf('     its = %i,  res = %i',its,res));
   if (res < err) 
      return
   end

   % fixed-point iteration solution update
   u = u + lam*f(1);
   v = v + lam*f(2);
   w = w + lam*f(3);

   % compute residuals
   f = chem_residual;
   res = max(abs(f));

end

   % declare the nested function chem_residual
   function resid = chem_residual
      resid(1) = k2*v + k1*u*(u+v+w-1);
      resid(2) = k1*u*(1-u-v-w) - (k2+k3)*v + k4*w;
      resid(3) = k3*v - k4*w;
   end

return
end
% end function
