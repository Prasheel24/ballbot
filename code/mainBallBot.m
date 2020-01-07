clc;
close all;
clear all;

% Initialise the symbols to be used
syms a;  %  alpha
syms b;  %  beta
syms g;  %  gamma
syms phi;  %  body angle
syms theta;  %  ball angle
syms phid;  %  body angular velocity
syms thetad;  %  ball angular velocity
syms Dc;  %  coloumb friction coefficient
syms Dv;  %  viscous friction coefficient
syms r;  %  ball radius
syms t;  %  Torque T

%Defining the Mass matrix
M = [a a+b*cosd(phi); a+b*cosd(phi) a+g+2*b*cosd(phi)];
%disp(M);

%Coriolis force Matrix
C = [-b*sind(phi)*phid*phid; -b*sind(phi)*phid*phid];
%disp(C);

%Potential energy Matrix
G = [0;(-b*g*sind(phi))/r];
%disp(G);

%Friction Matrix
D = [Dc*sign(thetad)+Dv*thetad; 0];
%disp(D);

%Torque Matrix
T = [t;0];

%The dynamic equation solving for theta'' and phi''
qdd = mtimes(inv(M),T) - mtimes(inv(M),G) - mtimes(inv(M),C); %-mtimes(inv(M),D)
disp(qdd);

%To define the A Matrix
dftheta = diff(qdd, theta);
%disp(dftheta);

dfphi = diff(qdd, phi);
%disp(dfphi);

dfthetad = diff(qdd, thetad);
%disp(dfthetad);

dfphid = diff(qdd, phid);
%disp(dfphid);

dfthetatd = diff(qdd, t);
%disp(dfthetatd);

%Calculate the A matrix with the parameters from the table in report
theta = 0;
phi = 0;
thetad = 0;
phid = 0;
% alpha = Iw + (mw+mb)*rw^2
a = 0.0174 + (2.44 + 51.66)*(0.1058)*(0.1058);
% beta = mb*rw*lb
b = 51.66 * 0.1058 * 0.69;
% gamma = Ib + mbIb^2
g = 12.48 + 51.66 * 0.69 * 0.69;
r = 0.1058;
N1 = (b * g * (g + b));
D1 = ( r * (b * b - a * g));
N2 = (a * b * g);
D2 = ( r * (-b * b + a * g));

N3 = -(a + g + 2 * b);
D3 = (b * b - a * g);
N4 = (a + b);
D4 = (b * b - a * g);

%To define B matrix
dfthetadd = N1/D1;
dfphidd = N2/D2;

dfthetatd1 = N3/D3;
dfphitd = N4/D4;

%Procedure to find A, B, C, D matrices to convert into State space
A = [0 0 1 0; 0 0 0 1; 0 dfthetadd 0 0; 0 dfphidd 0 0];
disp(A);
B = [0; 0; dfthetatd1; dfphitd];
disp(B);
C = eye(4);

%For the system to be controllable, the rank must be equal to matrix A rank
if length(A) == rank(ctrb(A, B))
    disp('The system is Controllable');
    %disp(rank(ctrb(A, B)));
else
    disp("System is not Controllable");
end
%For the system to be controllable, the rank must be equal to matrix A rank
if length(A) == rank(obsv(A, C))
    disp('The system is Observable');
    disp(obsv(A,C));
else
    disp("System is not Observable");
end

    
t0 = 0;
tf = 15;
tm = (t0+tf) / 2;
t = t0:0.1:tf; 
phiFunc = 0.0061*sech((9*(2*t-tm))/tm)-0.0061*sech(9*((2*t-15-tm)/(15-tm)));
plot(t, phiFunc)

% hold on
% 
% t0 = 0;
% tf = 19.54;
% tm = (t0+tf) / 2;
% t = t0:0.1:tf; 
% phiFuncFinal = 0.4062*sech((9*(2*t-tm))/tm)-0.4087*sech(9*((2*t-15-tm)/(15-tm)));
% plot(t, phiFuncFinal)
% hold off
grid on