%HW4.  Numerical PDEs APPM 6610.
%Pseudospectral method for wave eqn u_tt = (a(x)u_x)_x.%

N = 100;  %No. grid points in x.
% Establish correct indices of Fouurier transform for use later.
if (mod(N,2)==0)
    k = -N/2:N/2-1;
else
    k = -(N-1)/2:(N-1)/2;
end
tmax = 20;  %Maximum time of simulation.
Nt = 10000;  %No. of time points.
t = linspace(0,tmax,Nt);
h = t(2)-t(1);

%Solve ODE vj'' = (iax-a)vj for each value xj.
x = linspace(0,2*pi,N);
a = a_fun(x)';
ax = ax_fun(x)';
phi = phi_fun(x)';
psi = psi_fun(x)';

Uhat = zeros(N,Nt);  %Solution (in Fourier space) matrix.
U = zeros(N,Nt);  %Solution matrix.
cnzero = symfft(phi);   %Initial condition for ode u(x,0).
cnprimezero = symfft(psi);  %Initial condition for ode u_t(x,0).

%Solve ODE in time for the evolution of each grid point using RK4.
v = zeros(1,Nt);
for i = 1:N
   M = [0 1;1i*(k(i))*ax(i)-(k(i))^2*a(i) 0]; 
   y = [cnzero(i); cnprimezero(i)];
   v(1) = y(1);
   for j = 1:Nt-1
      k1 = M*y;
      k2 = M*(y+0.5*h*k1);
      k3 = M*(y+0.5*h*k2);
      k4 = M*(y+h*k3);
      y = y + (1/6)*h*(k1+2*k2+2*k3+k4);
      v(j+1) = y(1);
   end
   Uhat(i,:) = v;
end
for j = 1:Nt-1
    U(:,j) = symifft(Uhat(:,j));
end
%Plot solution.
figure()
%hold on
for i = 1:100:Nt
    plot(x,real(U(:,i)),'k','linewidth',1)
    hold on
    axis([0 7 -1 1])
    pause(0.1)
end
xlabel('x')

%Plot functions.
figure()
plot(x,a_fun(x),'k','linewidth',2)
xlabel('x')
legend('a(x)')
figure()
plot(x,phi_fun(x),'k','linewidth',2)
xlabel('x')
legend('phi(x)')
figure()
plot(x,psi_fun(x),'k','linewidth',2)
xlabel('x')
legend('psi(x)')

%%%%%%% Function Definitions %%%%%%%

%Symmetric fft.
function out = symfft(x)
    N = length(x);
    if (mod(N,2)==0)
        k = -N/2:N/2-1;
    else
        k = -(N-1)/2:(N-1)/2;
    end
    X = fft(x);
    out = fftshift(X);
end
%Symmetric inverse fft.
function out = symifft(x)
    N = length(x);
    if (mod(N,2)==0)
        k = -N/2:N/2-1;
    else
        k = -(N-1)/2:(N-1)/2;
    end
    X = ifftshift(x);
    out = ifft(X);
end
%Define function a(x).
function out = a_fun(x)
    out = ones(1,length(x));
end
%Define function a_x(x).
function out = ax_fun(x)
    out = zeros(1,length(x));
end
%Define function phi(x).
function out = phi_fun(x)
    out = exp(-10*(x-pi/2).^2);
    %out = cos(x);
end
%Define function psi(x).
function out = psi_fun(x)
    out = zeros(1,length(x));
end