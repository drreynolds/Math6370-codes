% Inputs the file u_sol.txt to determine problem/parallelism information.
% Inputs the files u_sol.000.000, u_sol.001.000, ..., to read the output from
% each MPI process at each time step.  Concatenates these solutions together
% for each test.
%
% plots the solution difference from the original (serial_f), and outputs
% the difference norm.

clear

% time step output to check
it = 5;

%%%% serial_f %%%%
cd serial_f;

% input general problem information
load u_sol.txt;
nx = u_sol(1);
ny = u_sol(2);
px = u_sol(3);
py = u_sol(4);
nt = u_sol(5);

% set local subdomain sizes (all but last in each dimension)
nloc = floor(nx/px);
mloc = floor(ny/py);

% allocate solution array
u = zeros(nx,ny);

% set time string
tstr = [ '.00' , sprintf('%i',it) ];
pname = [ 'sol.00' , sprintf('%i',it) , '.png' ];

% iterate over processors
for p=1:px*py
   
   % input this processes' output
   if (p < 11) 
      infile = [ 'u_sol.00' , sprintf('%i',p-1) , tstr ];
   elseif (p < 101)
      infile = [ 'u_sol.0' , sprintf('%i',p-1) , tstr ];
   else
      infile = [ 'u_sol.' , sprintf('%i',p-1) , tstr ];
   end
   data = load(infile);
   
   % extract this processes' subdomain information
   locN = data(1);
   locM = data(2);
   pcoords = data(3:4);
   t = data(5);
   
   % determine location of the processes' solution in overall domain
   is = nloc*pcoords(1) + 1;
   ie = is + locN - 1;
   js = mloc*pcoords(2) + 1;
   je = js + locM - 1;

   % read subdomain information into overall solution
   count = 6;
   for j=1:locM
      for i=1:locN
	 u(is+i-1,js+j-1) = data(count);
	 count = count + 1;
      end 
   end
   
   % clear subdomain information
   clear data
end

% store result as utrue;
utrue = u;
cd ..;



%%%% mpi_f %%%%
cd mpi_f;

% input general problem information
load u_sol.txt;
nx = u_sol(1);
ny = u_sol(2);
px = u_sol(3);
py = u_sol(4);
nt = u_sol(5);

% set local subdomain sizes (all but last in each dimension)
nloc = floor(nx/px);
mloc = floor(ny/py);

% clear solution array
u = zeros(nx,ny);

% iterate over processors
for p=1:px*py
   
   % input this processes' output
   if (p < 11) 
      infile = [ 'u_sol.00' , sprintf('%i',p-1) , tstr ];
   elseif (p < 101)
      infile = [ 'u_sol.0' , sprintf('%i',p-1) , tstr ];
   else
      infile = [ 'u_sol.' , sprintf('%i',p-1) , tstr ];
   end
   data = load(infile);
   
   % extract this processes' subdomain information
   locN = data(1);
   locM = data(2);
   pcoords = data(3:4);
   t = data(5);
   
   % determine location of the processes' solution in overall domain
   is = nloc*pcoords(1) + 1;
   ie = is + locN - 1;
   js = mloc*pcoords(2) + 1;
   je = js + locM - 1;
   
   % read subdomain information into overall solution
   count = 6;
   for j=1:locM
      for i=1:locN
	 u(is+i-1,js+j-1) = data(count);
	 count = count + 1;
      end 
   end
   
   % clear subdomain information
   clear data
end

% compute solution difference, norm
udiff = utrue - u;
udiff_norm = sqrt(sum(sum(abs(udiff.^2)))/nx/ny);
fprintf('mpi_f:  solution difference norm = %g, plotting.\n',udiff_norm);

% plot solution difference
xvals = linspace(0,1,nx);
yvals = linspace(0,1,ny);
h = surf(yvals,xvals,udiff);
shading flat
view([50 44])
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
title(sprintf('solution difference'),'FontSize',14)
cd ..;

disp('paused: hit enter to continue')
pause
   


%%%% mpi_c %%%%
cd mpi_c;

% input general problem information
load u_sol.txt;
nx = u_sol(1);
ny = u_sol(2);
px = u_sol(3);
py = u_sol(4);
nt = u_sol(5);

% set local subdomain sizes (all but last in each dimension)
nloc = floor(nx/px);
mloc = floor(ny/py);

% clear solution array
u = zeros(nx,ny);

% iterate over processors
for p=1:px*py
   
   % input this processes' output
   if (p < 11) 
      infile = [ 'u_sol.00' , sprintf('%i',p-1) , tstr ];
   elseif (p < 101)
      infile = [ 'u_sol.0' , sprintf('%i',p-1) , tstr ];
   else
      infile = [ 'u_sol.' , sprintf('%i',p-1) , tstr ];
   end
   data = load(infile);
   
   % extract this processes' subdomain information
   locN = data(1);
   locM = data(2);
   pcoords = data(3:4);
   t = data(5);
   
   % determine location of the processes' solution in overall domain
   is = nloc*pcoords(1) + 1;
   ie = is + locN - 1;
   js = mloc*pcoords(2) + 1;
   je = js + locM - 1;
   
   % read subdomain information into overall solution
   count = 6;
   for j=1:locM
      for i=1:locN
	 u(is+i-1,js+j-1) = data(count);
	 count = count + 1;
      end 
   end
   
   % clear subdomain information
   clear data
end

% compute solution difference, norm
udiff = utrue - u;
udiff_norm = sqrt(sum(sum(abs(udiff.^2)))/nx/ny);
fprintf('mpi_c:  solution difference norm = %g, plotting.\n',udiff_norm);

% plot solution difference
xvals = linspace(0,1,nx);
yvals = linspace(0,1,ny);
h = surf(yvals,xvals,udiff);
shading flat
view([50 44])
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
title(sprintf('solution difference'),'FontSize',14)
cd ..;

disp('  ')
disp('finished')
disp('  ')


% end of script


