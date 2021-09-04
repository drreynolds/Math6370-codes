% Daniel R. Reynolds
% SMU Mathematics
% Math 4370/6370
% 7 February 2015
%-----------------------------------------------------------------
% Description: 
%    Evolves the first-order 2D wave equations in time.
%=================================================================

% clear all old variables from memory
clear

stime = cputime;
% read problem parameters from input file (should be in this order):
%    nx - number of nodes in x-direction
%    ny - number of nodes in y-direction
%    nt - number of time steps
%    tstop - final time (will stop at nt or stop, whichever is 1st)
%    c - wave speed
%    dtoutput - time frequency for outputting solutions
vars = load('input_matlab.txt');
nx = vars(1);
ny = vars(2);
nt = vars(3);
tstop = vars(4);
c = vars(5);
dtoutput = vars(6);
ftime = cputime;
iotime = ftime-stime;

% declare solution variables on staggered grid:
%    u(i,j)  corresponds to u(x_i,y_j)
%    v1(i,j) corresponds to u_t(x_i,y_j)
%    v2(i,j) corresponds to c*u_x(x_{i-1/2},y_j)
%    v3(i,j) corresponds to c*u_y(x_i,y_{j-1/2})
stime = cputime;
u = zeros(nx,ny);
v1 = zeros(nx,ny);
v2 = zeros(nx,ny);
v3 = zeros(nx,ny);

% set grid spacing and mesh points
dx = 1/nx;
dy = 1/ny;
xspan_c = linspace(dx,1,nx);          % x-grid points at whole indices
xspan_h = linspace(dx/2,1-dx/2,nx);   % x-grid points at half indices
yspan_c = linspace(dy,1,ny);          % y-grid points at whole indices
yspan_h = linspace(dy/2,1-dy/2,ny);   % y-grid points at half indices


% set initial condition for solution and derivatives
for j=1:ny

   % y locations
   y_c = yspan_c(j);
   y_h = yspan_h(j);

   for i=1:nx

      % x locations
      x_c = xspan_c(i);
      x_h = xspan_h(i);
      
      % set initial conditions on u_x, u_y [gaussian blob]
      %    u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
      %    c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
      %    c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2))
      u(i,j) = exp(-100*((x_c-1/3)^2+(y_c-1/2)^2));
      v1(i,j) = 0;
      v2(i,j) = -200*c*(x_h-1/3)*exp(-100*((x_h-1/3)^2+(y_c-1/2)^2));
      v3(i,j) = -200*c*(y_h-1/2)*exp(-100*((x_c-1/3)^2+(y_h-1/2)^2));
      
   end 
end
ftime = cputime;
inittime = ftime-stime;


% set initial time, output initial solution
t = 0;
toutput = 0;
noutput = 0;
outname = 'u_sol.0000';
stime = cputime;
save(outname,'u');


% plot initial state 
figure(1), surf(yspan_c,xspan_c,u)
shading flat, view([50 44])
axis([0, 1, 0, 1, -1, 1])
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
title(sprintf('u(x,y) at t = %g,  mesh = %ix%i',t,nx,ny),'FontSize',16)
ftime = cputime;
iotime = iotime+ftime-stime;


% start time stepping
runtime = 0;
for it=1:nt

   % start timer
   stime = cputime;

   % set time step
   dt = min([dx/c/50, dy/c/50]);
   
   % compute solutions at new times
   for j=1:ny
      for i=1:nx

	 % first update v1 to get to half time step
	 %    get relevant values for this location
	 if (i==nx)
	    v2_E = v2(1,j);
	 else
	    v2_E = v2(i+1,j);
	 end
	 v2_W = v2(i,j);
	 if (j==ny) 
	    v3_N = v3(i,1);
	 else
	    v3_N = v3(i,j+1);
	 end
	 v3_S = v3(i,j);
	 %    update v1
	 v1(i,j) = v1(i,j) + c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);
	 
      end
   end

   for j=1:ny
      for i=1:nx
	 % next update v2 & v3 to get to full time step
	 %    get relevant values for this location
	 if (i==1)
	    v1_W = v1(nx,j);
	 else
	    v1_W = v1(i-1,j);
	 end
	 v1_E = v1(i,j);
	 if (j==1)
	    v1_S = v1(i,ny);
	 else
	    v1_S = v1(i,j-1);
	 end
	 v1_N = v1(i,j);
	 % update v2 and v3
	 v2(i,j) = v2(i,j) + c*dt/dx*(v1_E - v1_W);
	 v3(i,j) = v3(i,j) + c*dt/dy*(v1_N - v1_S);
      
      end
   end
   
   % update and plot solution
   u = u + dt*v1;

   % update time
   t = t + dt;
%   disp(sprintf(' time step %i: dt = %g, t = %g',it,dt,t));
   ftime = cputime;
   runtime = runtime+ftime-stime;

   % stop simulation if we've reached tstop
   if (t >= tstop) 
      break;
   end
   
   % output solution periodically
   if (abs(t-toutput-dtoutput) <= 100*eps)
      stime = cputime;
      toutput = t;
      noutput = noutput+1;
      if (noutput < 10)
	 outname = [ 'u_sol.000', num2str(noutput) ];
      elseif (noutput < 100)
	 outname = [ 'u_sol.00', num2str(noutput) ];
      elseif (noutput < 1000)
	 outname = [ 'u_sol.0', num2str(noutput) ];
      else
	 outname = [ 'u_sol.', num2str(noutput) ];
      end
      save(outname,'u');
   
      figure(1), surf(yspan_c,xspan_c,u)
      shading flat, view([50 44])
      axis([0, 1, 0, 1, -1, 1])
      xlabel('x','FontSize',14), ylabel('y','FontSize',14)
      title(sprintf('u(x,y) at t = %g,  mesh = %ix%i',t,nx,ny),'FontSize',16)
      ftime = cputime;
      iotime = iotime+ftime-stime;

   end

end


% output final solution 
stime = cputime;
toutput = t;
noutput = noutput+1;
if (noutput < 10)
   outname = [ 'u_sol.000', num2str(noutput) ];
elseif (noutput < 100)
   outname = [ 'u_sol.00', num2str(noutput) ];
elseif (noutput < 1000)
   outname = [ 'u_sol.0', num2str(noutput) ];
else
   outname = [ 'u_sol.', num2str(noutput) ];
end
save(outname,'u');

figure(1), surf(yspan_c,xspan_c,u)
shading flat, view([50 44])
axis([0, 1, 0, 1, -1, 1])
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
title(sprintf('u(x,y) at t = %g,  mesh = %ix%i',t,nx,ny),'FontSize',16)
ftime = cputime;
iotime = iotime+ftime-stime;


% output runtime
disp(sprintf(' total initialization time = %g', inittime));
disp(sprintf(' total input/output time   = %g', iotime));
disp(sprintf(' total simulation time     = %g', runtime));


% end program
