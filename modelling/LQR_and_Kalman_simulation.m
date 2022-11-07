%{
    Author: Andrea Nicetto
    Project: Maglev2022 - https://github.com/Nicetto/Maglev2022
%}
clc;clear;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");

if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

%% Linearization of system

x_start = zeros(12,1); x_start(3) = eq;
sys = maglevSystem(x_start, params, approximationType);

uLp = zeros(params.solenoids.N,1); % inputs
xLp = zeros(12,1); xLp(3) = eq; % equilibrium point around which to linearize the system

delta = 1e-4;
dimX = 12;
dimU = params.solenoids.N;
dimY = 3*length(params.sensor.x); 

h = waitbar(0);
max_dim = dimX+dimU+dimX;

% d(x/dx) empty matrix of dimension x*x (12x12) where
% on rows there's df(i=1:12)/.
% on columns there's ./dx(i=1:12)

A = zeros(dimX,dimX);
for i = 1:dimX % for every column
    % linearize the entire column
    A(:,i) = (sys.f(xLp+delta*(i==1:dimX)',uLp) ...
             -sys.f(xLp-delta*(i==1:dimX)',uLp)) ...
             /(2*delta);
    waitbar(i/max_dim);
end

% d(x/du) empty matrix of dimension x*u (12x4)
B = zeros(dimX,dimU);  
for i = 1:dimU
    B(:,i) = (sys.f(xLp,uLp+delta*(i==1:dimU)') ...
             -sys.f(xLp,uLp-delta*(i==1:dimU)')) ...
             /(2*delta);
    waitbar((i+dimX)/max_dim);
end

% d(y/dx) empty matrix of dimension y*x (9x12)
C = zeros(dimY,dimX);
for i = 1:dimX
    C(:,i) = (sys.h(xLp+delta*(i==1:dimX)', uLp) ...
             -sys.h(xLp-delta*(i==1:dimX)', uLp)) ...
             /(2*delta);
    waitbar((i+dimX+dimU)/max_dim);
end
close(h);

% d(y/du) = 0
D = zeros(dimY, dimU);

%% LQR controller

% - Two of the states are uncontrollable, so have to disregard those
Ared = A([1,2,3,4,5,7,8,9,10,11],[1,2,3,4,5,7,8,9,10,11]);
Bred = B([1,2,3,4,5,7,8,9,10,11],:);
Cred = C(:,[1,2,3,4,5,7,8,9,10,11]);
Dred = zeros(15,4);

%Setting the parameters for the LQR controller
Q = diag([1e5,1e5,1e2, 1e1,1e1, 1e2,1e2,1e1, 1e2,1e2]);
R = 1e-0*eye(4);

%LQR controller 
[Kred, S, poles] = lqr(Ared,Bred,Q,R);

%Adding back the zero columns
K = [Kred(:,1:5), zeros(4,1), Kred(:,6:end), zeros(4,1)];

%Simulation of the system with eh LQR controller 
% starting from the initial position x0 
x0 = xLp+[0.0006,0.0002,0.0008, -0.0002,0,0, 0,0,0, 0,0,0]';
tspan  = linspace(0,3,100);

odefun = @(t,x) sys.f(x,-K*(x-xLp)-uLp);
[t,x] = ode45(odefun, tspan, x0);

%% Kalman filter

%Defining covariances
Q_kalman = 10*eye(4);

R_kalman = 0.05*eye(15);

%Implementing a Kalman filter
[kalmf,L,P] = kalman(ss(Ared,Bred,Cred,[]),Q_kalman,R_kalman);

%Get the measurements from the sensors
n = length(tspan);

yLp=sys.h(xLp,zeros(4,1)); %Remember that the Kalman filter is applied around xLp

y = zeros(n,15);
for i = 1:size(x,1)
    y(i,:) = sys.h(x(i,:)',-K*(x(i,:)'-xLp));
end

%Estimating the position of the levitating magnet with the Kalman filter
uInterp = @(t) interp1(tspan, (-K*(x'-xLp))', t)';
yInterp = @(t) interp1(tspan, y, t)';

odefun = @(t,xhat) Ared*xhat + Bred*uInterp(t) + L*(yInterp(t)-yLp - Cred*(xhat-xLp([1:5,7:11])));

x0red = x0([1:5,7:11]);
[t,xhat] = ode45(odefun,tspan,x0red);

%Compute the error between the estimated trajectory and the real one
for k = 1:100
   estimation_error(k,1)=norm(xhat(k,1:3)-x(k,1:3));
end

%% Plotting  results

% z evolution
figure();
clf; grid on; hold on;
xlabel('$t$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('z', 'fontsize', 16, 'interpreter', 'latex')
plot(t,x(:,3), 'b', 'linewidth', 2);

% y evolution
figure();
clf; grid on; hold on;
xlabel('$t$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('y', 'fontsize', 16, 'interpreter', 'latex')
plot(t,x(:,2), 'r', 'linewidth', 2);

% x evolution
figure();
clf; grid on; hold on;
xlabel('$t$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('x', 'fontsize', 16, 'interpreter', 'latex')
plot(t,x(:,1), 'r', 'linewidth', 2);

% Other coordinates
figure();
clf; grid on; hold on;
xlabel('$t$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$\eta$', 'fontsize', 16, 'interpreter', 'latex')
plot(t,x(:,4:6), 'r', 'linewidth', 2);

% plotting the estimated error
figure();
clf; grid on; hold on;
xlabel('$t$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$||x-\hat{x}||$', 'fontsize', 16, 'interpreter', 'latex')
plot(t(2:50,1),estimation_error(2:50,1),'r-','linewidth',2)


