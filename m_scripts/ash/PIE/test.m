%{
% TyProg020101.m - Solve KdV eq. u_t + uu_x + u_xxx = 0 on [-pi,pi] by 
% FFT WITHOUT integrating factor used by Trefethen in p27.m.
% The effect of this is the requirement for a greatly reduced 
% time step to retain stability.
% Set up grid and two-soliton initial data:

N = 256;
x = (2*pi/N)*(-N/2:N/2-1)';
A = 25; B = 16; clf, drawnow, set(gcf,'renderer','zbuffer')
u = 3*A^2*sech(.5*(A*(x+2))).^2 + 3*B^2*sech(.5*(B*(x+1))).^2; 
v = fft(u); k = [0:N/2-1 0 -N/2+1:-1]'; ik3 = 1i*k.^3;

% The time step. Original from p27.m was dt = .4/N^2;
dt = 1.3864e-06; % Good
%dt = 1.3865e-06; % Start to see noise at end of time integration 
%dt = 1.3866e-06; % Bad

% Solve PDE and plot results:
tmax = 0.006;
nplt = floor((tmax/25)/dt); 
nmax = round(tmax/dt); 
udata = u; 
tdata = 0; 

h = waitbar(0,'please wait...'); 
for n = 1:nmax
    t = n*dt; 
    g = -.5i*dt*k; g2 = -2*k.^2; 
    a = g.*(fft((ifft(v)).^2) + g2.*v);
    b = g.*(fft((ifft(v+a/2)).^2) + g2.*(v+a/2)); 
    c = g.*(fft((ifft(v+b/2)).^2) + g2.*(v+b/2)); 
    d = g.*(fft((ifft(v+c)).^2) + g2.*(v+c));
    v = v + (a + 2*(b+c) + d)/6;
    if mod(n,nplt) == 0
        u = real(ifft(v)); waitbar(n/nmax)
        udata = [udata u]; tdata = [tdata t]; 
    end
    % 4th-order
    % Runge-Kutta
end

waterfall(x,tdata,udata'), colormap(1e-6*[1 1 1]); 
view(-20,25)
xlabel x, ylabel t, axis([-pi pi 0 tmax 0 2000]), grid, 'off', set(gca,'ztick',[0 2000]), close(h), pbaspect([1 1 .13]);
%}


%{
% Approximates the derivative of a periodic function f(x) using Fourier
% transforms.  Requires a linear discretization and a domain (a,b] such
% that f(x) = f(x+b-a)
 
% domain
a = 1;
b = 1 + pi/2;
N = 100;
dx = (b-a)/N;
x = a + dx*(0:N-1);
 
% function
w = 2;
f = sin(w*x).^2;
 
% exact derivatives
dfdx = 2*w*sin(w*x).*cos(w*x);
d2fdx2 = 4*w^2*cos(w*x).^2 - 2*w^2;
 
% fourier derivatives
%Nx = size(x,2);
%k = 2*pi/(b-a)*[0:Nx/2-1 0 -Nx/2+1:-1];
%dFdx = ifft(1i*k.*fft(f));
%d2Fdx2 = ifft(-k.^2.*fft(f));

%modified SOH
Nx = length(x);
dk = 2*pi/(Nx*dx);
k= fftshift(dk*[-Nx/2 : Nx/2-1]);
dFdx = ifft(1i*k.*fft(f));
d2Fdx2 = ifft(-k.^2.*fft(f));

% graph result
clf;
plot(x,dfdx,'r-',x,d2fdx2,'g-',x,dFdx,'k:',x,d2Fdx2,'b:','LineWidth',2);
legend('df/dx','d^2f/dx^2','Fourier df/dx','Fourier d^2f/dx^2');
%}

%{
%Solving Heat Equation using pseudo-spectral and Forward Euler
%u_t= \alpha*u_xx
%BC= u(0)=0, u(2*pi)=0
%IC=sin(x)

%Grid
N = 64; %Number of steps
h = 2*pi/N; %step size
x = h*(1:N); %discretize x-direction

alpha = .5; %Thermal Diffusivity constant
t = 0;
dt = .001;

%Initial conditions
v = sin(x);
v_o=v;
%k=(1i*[0:N/2-1 0 -N/2+1:-1]); %original
k=fftshift([-32:31]);  %my version
ik2=(1i*k).^2;

%Setting up Plot
tmax = 5; tplot = .1;
plotgap= round(tplot/dt);
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; 
tdata = t;

cnt = 0;
v_hat = fft(v); %Fourier Space

for i = 1:nplots
    for n = 1:plotgap
        v_hat = v_hat+dt*alpha*ik2.*v_hat; %FE timestepping
        cnt = cnt+1;
    end
    v = real(ifft(v_hat)); %Back to real space
    data(i+1,:) = v;
    t=t+plotgap*dt;
    tdata = [tdata; t]; %Time vector
end

v_final=v;

%Plot using mesh
mesh(x,tdata,data), grid on,
view(-60,55), xlabel x, ylabel t, zlabel u, zlabel u

%SOH trials
v_hat = fft(v_o); %Fourier Space
for i = 1:(tmax/dt)
    v_hat = v_hat + dt*alpha*ik2.*v_hat; %FE timestepping
end
v_soh = real(ifft(v_hat)); %Back to real space
t=i*dt;

hold on; plot3(x, ones(size(x))*t, v_final, 'k', 'linewidth', 1); hold off;
%}

%{
%Solving Heat Equation using pseudospectral methods with Backwards Euler:
%u_t= \alpha*u_xx
%BC = u(0)=0 and u(2*pi)=0 (Periodic)
%IC=sin(x)
clear all; clc;
 
%Grid
N = 64; h = 2*pi/N; x = h*(1:N);
 
% Initial conditions
v = sin(x);
alpha = .5;
t = 0;
dt = .001; %Timestep size
 
%(ik)^2 Vector
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k=1i*fftshift([-32:31]);  %my version

k2=k.^2;
 
%Setting up Plot
tmax = 5; tplot = .1;
plotgap= round(tplot/dt);
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;
 
 
for i = 1:nplots
    v_hat = fft(v); %Converts to fourier space
    for n = 1:plotgap
        v_hat = v_hat./(1-dt*alpha*k2); %Backwards Euler timestepping
    end
    v = ifft(v_hat); %Converts back to real Space
    data(i+1,:) = real(v); %Records data
    t=t+plotgap*dt; %Records time
    tdata = [tdata; t];
end
 
%Plot using mesh
mesh(x,tdata,data), grid on, %axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel u, zlabel u,
%}


%{
%Solving 1D Allen-Cahn Eq using pseudo-spectral and Implicit/Explicit method
%u_t=u_{xx} + u - u^3
%where u-u^3 is treated explicitly and u_{xx} is treated implicitly
%BC = u(0)=0, u(2*pi)=0 (Periodic)
%IC=.25*sin(x);
clear all; clc;

%Grid and Initial Data
N = 8000; h = 2*pi/N; x = h*(1:N); t = 0;

dt = .001; %timestep size
epsilon= .001;

%initial conditions
v = .25*sin(x);

%(ik) and (ik)^2 vectors
k=(1i*[0:N/2-1 0 -N/2+1:-1]);
k2=k.^2;

%setting up plot
tmax = 5; tplot = .2;
plotgap= round(tplot/dt);
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; tdata = t;
  
for i = 1:nplots
    for n = 1:plotgap
        v_hat = fft(v); %converts to Fourier space
        vv = v.^3; %computes nonlinear term in real space
        vv = fft(vv); %converts nonlinear term to Fourier space
        v_hat = (v_hat*(1/dt+1) - vv)./(1/dt-k2*epsilon); %Implicit/Explicit
        v = ifft(v_hat); %Solution back to real space
    end
    data(i+1,:) = real(v); %Records data each "plotgap"
    t=t+plotgap*dt; %Real time
    tdata = [tdata; t];
end

mesh(x,tdata,data), grid on, axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel u
%}


%%{
%Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
%u_t= epsilon(u_{xx}+u_{yy}) + u - u^3
%where u-u^3 is treated explicitly and epsilon(u_{xx} + u_{yy}) is treated implicitly
%BC = Periodic
%IC=v=sin(2*pi*x)+0.001*cos(16*pi*x;
clear all; clc;

%Grid
N = 256; h = 1/N; x = h*(1:N);
dt = .01;

%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);

%initial conditions
%v=sin(2*pi*xx)+0.001*cos(16*pi*xx);
v = rand(size(xx))*2-1;
epsilon=.01;

%(ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;

[kxx,kyy]=meshgrid(k2x,k2y);
        
for n = 1:500
    v_nl=v.^3; %calculates nonlinear term in real space
    %FFT for linear and nonlinear term
    v_nl = fft2(v_nl);
    v_hat=fft2(v);
    vnew=(v_hat*(1+1/dt)-v_nl)./ ...
       (-(kxx+kyy)*epsilon+1/dt); %Implicit/Explicit timestepping
    %converts to real space in x-direction
    v=ifft2(vnew);
    %Plots each timestep
    mesh(v); title(['Time ',num2str(n)]); axis([0 N 0 N -1 1]);
    xlabel x; ylabel y; zlabel u;
    view(43,22); drawnow;
end
%}