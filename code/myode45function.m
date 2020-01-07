% % function dynFuncTheta = myode45function(t,theta)
% % a = 0.0174 + (2.44 + 51.66)*(0.1058)*(0.1058);
% % b = 51.66 * 0.1058 * 0.69;
% % g = 12.48 + 51.66 * 0.69 * 0.69;
% % r = 0.1058;
% % 
% % t0 = 0;
% % tf = 15;
% % tm = (t0+tf) / 2;
% % phiFunc = 0.0061*sech((9*(2*t-tm))/tm)-0.0061*sech(9*((2*t-15-tm)/(15-tm)));
% % dynFuncTheta = [theta(2), -((a+g+2*b*cosd(phiFunc))/(a+b*cosd(phiFunc))) ];
% 
% syms a;
% syms b;
% syms g;
% syms phi;
% syms theta;
% syms phid;
% syms thetad;
% syms Dc;
% syms Dv;
% syms r;
% syms g;
% syms t;
% M = [a a+b*cosd(phi);a+b*cosd(phi) a+g+2*b*cosd(phi)];
% %disp(M);
% C = [-b*sind(phi)*phid*phid; -b*sind(phi)*phid*phid];
% %disp(C);
% G = [0;(-b*g*sind(phi))/r];
% %disp(G);
% D = [Dc*sign(thetad)+Dv*thetad; 0];
% %disp(D);
% 
% T = [t;0];
% qdd = mtimes(inv(M),T) - mtimes(inv(M),G) - mtimes(inv(M),C);
% %disp(qdd(1));
% thetadd = qdd(1);

% theta = 0;
% phi = 0;
% thetad = 0;
% phid = 0;
% t0 = 0;
% tf = 15;
% tm = (t0+tf) / 2;
% t = t010.011tf; 
% phiFunc = 0.0061*sech((9*(2*t-tm))/tm)-0.0061*sech(9*((2*t-15-tm)/(15-tm)));
% disp(phiFunc(1));
% plot(t, phiFunc)
% grid on



% [t,y] = RK45(@(t,y)[y(2); (1-y(1)^2)*y(2)-y(1)],[0 20],[2; 0]);
% subplot(1,2,2);
% plot(t,y,'-o')
% title('Solution of van der Pol Equation (\mu = 1) with RK45');
% xlabel('Time t');
% ylabel('Solution y');
% % legend('y_1','y_2');
function plot_vector_ode 
close all
theta0 = 0.075;
phi0 = 0.0975;x
tstart = 0;
tstop = 15;
theta0 = theta0/180;
phi0 = phi0 / 180;
%-((27*(a+b*cos(w(3))))/(a*g-b*b*cos(w(3))*cos(w(3)))) -((b*w(4)*2*sin(w(3))*(a+b*cos(w(3))))/(a*g-b*b*cos(w(3))*cos(w(3)))) +((b*w(4)*2*sin(w(3)))/(a*g-b*b*cos(w(3))*cos(w(3)))) +((a*b*g*sin(w(3)))/(r*(a*g-b*b*cos(w(3))*cos(w(3)))))
a = 0.0174 + (2.44 + 51.66)*(0.1058)*(0.1058);
b = 51.66 * 0.1058 * 0.69;
g = 12.48 + 51.66 * 0.69 * 0.69;
r = 0.1058;
t0 = 0;
tf = 15;
tm = (tf+t0)/2;
k = 9;
a1 = 0.0061;
a2 = -0.0061;

% fun = -(((a+g+2*b*cos(phi))*diff(phi))/(a+b*cos(phi))) + ((b*sin(phi)*((phi)^2))/(a+b*cos(phi))) + ((b*g*sin(phi))/(r*(a+b*cos(phi))));
w0 =  [theta0; 0; phi0; 0];
[times, sols] = ode45(@odefunc, [tstart,tstop], w0);
plot(sols(:,1),times);
disp(sols);
function dwdt = odefunc(t,w)
  theta = w(1);
  thetad = w(2);
  phi = w(3);
  phid = w(4);
  phidd = ((-a1*((2*k)/(tm-t0))^2) * (sech(k*(2*t-tm-t0)/tm-t0) - (2*((tanh(k*(2*t-tm-t0)/tm-t0))^2)*(sech(k*(2*t-tm-t0)/tm-t0)))) - (-a2*((2*k)/(tf-tm))^2) * (sech(k*(2*t-tf-tm)/tf-tm) - (2*((tanh(k*(2*t-tf-tm)/tf-tm))^2)*(sech(k*(2*t-tf-tm)/tf-tm)))));%-((20*(a+b*cos(w(3))))/((a*g)-(b*b*cos(w(3))*cos(w(3))))) -((b*w(4)*2*sin(w(3))*(a+b*cos(w(3))))/(a*g-b*b*cos(w(3))*cos(w(3)))) +((a*b*w(4)*2*sin(w(3)))/(a*g-b*b*cos(w(3))*cos(w(3)))) +((a*b*g*sin(w(3)))/(r*(a*g-b*b*cos(w(3))*cos(w(3)))));
  dwdt = [w(2);-(((a+g+2*b*cos(w(3)))*phidd)/(a+b*cos(w(3)))) + ((b*sin(w(3))*((w(3))^2))/(a+b*cos(w(3)))) + ((b*g*sin(w(3)))/(r*(a+b*cos(w(3))))); w(4); phidd];
end
end
% k = 9;
% tm = (tf+t0)/2;
% options = optimset('PlotFcns',@optimplotfval);
% t0 = 0;
% tf = 15;
% a1 = 0.0061;
% a2 = -0.0061;
 %fun = @(t)((-a1*((2*k)/(tm-t0))^2) * (sech(k*(2*t-tm-t0)/tm-t0) - (2*((tanh(k*(2*t-tm-t0)/tm-t0))^2)*(sech(k*(2*t-tm-t0)/tm-t0)))) - (-a2*((2*k)/(tf-tm))^2) * (sech(k*(2*t-tf-tm)/tf-tm) - (2*((tanh(k*(2*t-tf-tm)/tf-tm))^2)*(sech(k*(2*t-tf-tm)/tf-tm)))));% - (-a2*((2*k)/(tf-tm)).^2)*((sech((k*(2*t-tf-tm))/(tf-tm))) -2*((tanh((k*(2*t-tf-tm))/(tf-tm)))).^2)*(sech((k*(2*t-tf-tm))/(tf-tm))) );
% 
% 
% x = fminsearch(fun, [0, 15], options)
% 
% fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
% x0 = [-1.2,1];
% x = fminsearch(fun,x0,options)

%disp(fun);             
  %disp( );



