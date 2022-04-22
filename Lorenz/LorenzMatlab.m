clc
clear
close all

% Tarea
h = 1e-2;
tfin = 200;
N = ceil((tfin-h)/h);
t = h + (0:N)*h;

% Ecuaciones
% \dot{x}= sigma*(y-x)
% \dot{y}= x*(rho-z)-y
% \dot{z}= x*y -beta*z

% Parámetros
sigma = 10;
beta = 8/3;
rho = 28;

% Condiciones iniciales
x0 = 0.1;
y0 = 0.1;
z0 = 0.1;

%================== EULER ===================================

% prealojamos x, y & z
xe = [x0 zeros(1,N-1)];
ye = [y0 zeros(1,N-1)];
ze = [z0 zeros(1,N-1)];

tic
for n = 1:N
    xe(n+1) = xe(n) + h*(sigma * (ye(n) - xe(n)) );
    ye(n+1) = ye(n) + h*(xe(n) * (rho - ze(n)) - ye(n) );
    ze(n+1) = ze(n) + h*(xe(n) * ye(n) - (beta * ze(n)) );
end
toc

% graficamos sólo este método
% comentamos, ya que funciona
% para x
%figure
%plot(xe,ye)
%figure
%plot(xe,ze)
%figure
%plot(t,xe)

% para y
%figure
%plot(ye,xe)
%figure
%plot(ye,ze)
%figure
%plot(t,ye)

% para z
%figure
%plot(ye,xe)
%figure
%plot(ye,ze)
%figure
%plot(t,ye)

% todas juntas
%figure
%plot3(xe,ye,ze)
%%%%axis equal
%xlabel('x(t)')
%ylabel('y(t)')
%zlabel('z(t)')

%============== fin EULER ===================================

%================= RK4 ======================================
% prealojamos x & y
xRK4 = [x0 zeros(1,N-1)];
yRK4 = [y0 zeros(1,N-1)];
zRK4 = [z0 zeros(1,N-1)];

tic
for n = 1:N
    % \dot{x}= sigma*(y-x)
    k11 = sigma*(yRK4(n)-xRK4(n));
    k21 = sigma*(yRK4(n)-(xRK4(n) + 0.5*h*k11));
    k31 = sigma*(yRK4(n)-(xRK4(n) + 0.5*h*k21));
    k41 = sigma*(yRK4(n)-(xRK4(n) + h*k31));
    
    xRK4(n+1) =  xRK4(n) + (h/6)*(k11 + 2*k21 + 2*k31 + k41);
       
    % \dot{y}= x*(rho-z)-y
    k12 = xRK4(n)*(rho-zRK4(n))-yRK4(n);
    k22 = xRK4(n)*(rho-zRK4(n))-(yRK4(n) + 0.5*h*k12);
    k32 = xRK4(n)*(rho-zRK4(n))-(yRK4(n) + 0.5*h*k22);
    k42 = xRK4(n)*(rho-zRK4(n))-(yRK4(n) + h*k32);
    
    yRK4(n+1) = yRK4(n) + (h/6) * (k12 + 2*k22 + 2*k32 + k42);
    
    % \dot{z}= x*y -beta*z
    k13 = xRK4(n)*yRK4(n) - beta*zRK4(n);
    k23 = xRK4(n)*yRK4(n) - beta*(zRK4(n) + 0.5*h*k13);
    k33 = xRK4(n)*yRK4(n) - beta*(zRK4(n) + 0.5*h*k23);
    k43 = xRK4(n)*yRK4(n) - beta*(zRK4(n) + h*k33);
    
    zRK4(n+1) = zRK4(n) + (h/6) * (k13 + 2*k23 + 2*k33 + k43);   
end
toc

% graficamos sólo este método
% comentamos, ya que funciona
% para x
%figure
%plot(xRK4,yRK4)
%figure
%plot(xRK4,zRK4)
%figure
%plot(t,xRK4)

% para y
%figure
%plot(yRK4,xRK4)
%figure
%plot(yRK4,zRK4)
%figure
%plot(t,yRK4)

% para z
%figure
%plot(zRK4,xRK4)
%figure
%plot(zRK4,yRK4)
%figure
%plot(t,zRK4)

% todas juntas
%figure
%plot3(xRK4,yRK4,zRK4)
%%%%axis equal
%xlabel('x(t)')
%ylabel('y(t)')
%zlabel('z(t)')

%================= fin RK4 ==================================

%============= EULER vs RK4 =================================
% x,y,z (pares)
figure
tiledlayout(1,3)
ax1 = nexttile;
plot(xe,ye)
hold on
plot(xRK4,yRK4)
xlabel('x(t)')
ylabel('y(t)')
legend('Euler', 'RK4')
title(ax1, 'x,y')

ax2 = nexttile;
plot(xe,ze)
hold on
plot(xRK4,zRK4)
xlabel('x(t)')
ylabel('z(t)')
legend('Euler', 'RK4')
title(ax2, 'x,z')

ax3 = nexttile;
plot(ye,ze)
hold on
plot(yRK4,zRK4)
xlabel('y(t)')
ylabel('z(t)')
legend('Euler', 'RK4')
title(ax3, 'y,z')

% Temporales
figure
tiledlayout(1,3)
ax1 = nexttile;
plot(t,xe)
hold on
plot(t,xRK4)
xlabel('t')
ylabel('x')
legend('Euler', 'RK4')
title(ax1, 'x,y')

ax2 = nexttile;
plot(t,ye)
hold on
plot(t,yRK4)
xlabel('t')
ylabel('y')
legend('Euler', 'RK4')
title(ax2, 'x,z')

ax3 = nexttile;
plot(t,ze)
hold on
plot(t,zRK4)
xlabel('t')
ylabel('z')
legend('Euler', 'RK4')
title(ax3, 'y,z')

% Todas juntas
figure
plot3(xe,ye,ze)
hold on
plot3(xRK4,yRK4,zRK4)
%%%%axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
legend('Euler', 'RK4')



dt = 0.01;
%bx = benchmark(xe,xRK4,t, dt)
%by = benchmark(ye,yRK4,t, dt)
%bz = benchmark(ze,zRK4,t, dt)




