% Daniel R. Reynolds
% SMU Mathematics
% Math 4370/6370
% 7 February 2015
%-----------------------------------------------------------------
% Description: 
%    Computes pi through numerical integration via
%        pi = 4*int_0^1 1/(1+x^2) dx
%    We use a simple midpoint rule for integration, over 
%    subintervals of fixed size 1/n, where n is a user-input 
%    parameter.
%=================================================================

% clear all old variables from memory
clear

% declare variables to be used
n = 0;
i = 0;
h = 0;
x = 0;
f = 0;
a = 0;
stime = 0;
ftime = 0;
pi_true = pi;
pi = 0;


%  set the integrand function handle
f = @(a) 4.d0 / (1.d0 + a*a);

% input the number of intervals
n = input('Enter the number of intervals (0 quits): ');
if (n < 1) 
   return
end

%  start timer
stime = cputime;

% set subinterval width
h = 1/n;

% perform integration over n intervals
pi = 0;
for i=1:n
   x = h*(i - 0.5);
   pi = pi + h*f(x);
end

% stop timer
ftime = cputime;

% output computed value and error
disp(sprintf(' computed pi = %.12e',pi));
disp(sprintf('     true pi = %.12e',pi_true));
disp(sprintf('       error = %g',pi_true - pi));
disp(sprintf('     runtime = %.5e',ftime-stime));
  

% end program 
