%HW4.  Numerical PDEs APPM 6610.
%Finite Difference Method%

N = 99;  %No. grid points in x.
x = linspace(1/N,2*pi,N)';
tmax = 50;  %Maximum time of simulation.
Nt = 10000;  %No. of time points.
t = linspace(0,tmax,Nt);
hx = x(2)-x(1);
ht = t(2)-t(1);
cols = int16(tmax/ht);  %Number of time values for simulation.
U = zeros(N+1,cols);

%Evaluate functions on spacial grids.
phi = phi_fun(x);
psi = psi_fun(x);
U(2:end,1) = phi;
U(1,1) = U(end,1);  %periodic bc.
col = 1;
i = 1;
%initial condition.
u0 = zeros(N+1,1);  u1 = zeros(N+1,1);  u = zeros(N+1,1);

%Compute u(ht) = u1%.

%Form matrix A.
main_diag = zeros(1,N);
sup_diag = zeros(1,N-1);
for j = 1:N-1
    main_diag(j) = -1*(a_fun(x(j)-pi/N)+a_fun(x(j)+pi/N));
end
main_diag(N) = main_diag(1);
for j = 1:N-1
   sup_diag(j) = a_fun(x(j+1));
end
sub_diag = sup_diag;
A = diag(main_diag,0)+diag(sup_diag,1)+diag(sub_diag,-1);
A(1,N) = a_fun(pi/N);
A(N,1) = A(1,N);
B = eye(N) + 0.5*(ht/hx)^2*A;
u1(2:end) = B*(phi+ht*psi);
u1(1) = u1(end);
U(:,2) = u1;
u0 = u1;
col = 2;
i = 2;

%Explicit central difference scheme.
while (col < cols)
    t = i*ht;
    col = col + 1;
    %compute solution at time ht.
    u(2:end) = 2*B*u1(2:end)-u0(2:end);
    u(1) = u(end);
    U(:,col) = u;
    %Prepare for next loop.
    u0 = u1;    u1 = u;
    i = i + 1;
end

%animation
for i = 1:100:cols
    plot([0;x],U(:,i),'k','linewidth',1)
    hold on
    axis([0 7 -1 1])
    pause(0.2)
end
xlabel('x')

figure()
plot([0;x],phi_fun([0;x]),'k')

%%%%%%% Function Definitions %%%%%%%
%Define function a(x).
function out = a_fun(x)
    %out = 1;
    if (x<pi)
        out = 0.95;
    else
        out = 1;
    end
    
end
%Define function phi(x).
function out = phi_fun(x)
    out = exp(-10*(x-pi/2).^2);
    %out = cos(x);
end
%Define function psi(x).
function out = psi_fun(x)
    out = zeros(length(x),1);
end







