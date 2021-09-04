% Daniel R. Reynolds
% SMU, Mathematics
% Math 4370/6370
% 7 February 2015
%-----------------------------------------------------------------
% Description: 
%    Computes the dot-product of two vectors.
%=================================================================

% clear all old variables from memory
clear

% declare variables to be used
n = 100000;
i = 0;
a = zeros(n,1);
b = zeros(n,1);
sum = 0;
stime = 0;
ftime = 0;

% initialize a and b
for i=1:n
  a(i) = 0.001*i/n;
  b(i) = 0.001*(n-i)/n;
end

% start timer
stime = cputime;

% compute dot-product
sum = 0.0;
for i=1:n
   sum = sum + a(i)*b(i);
end

% stop timer
ftime = cputime;

% output computed value and runtime
disp(sprintf(' dot-product = %g',sum));
disp(sprintf('     runtime = %g',ftime-stime));

% end program
