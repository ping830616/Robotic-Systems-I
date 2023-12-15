clear all
syms th1 th1d th2 th2d 
delete("ode1.m")
delete("ode1_block.slx")
g=9.81; 
m1=0.5;
m2=0.5;
a1=2;
a2=2;
L=(m2*((a1*cos(th1)*th1d + a2*cos(th1)*cos(th2)*th1d - a2*sin(th1)*sin(th2)*th2d)^2 + (a1*sin(th1)*th1d + a2*cos(th2)*sin(th1)*th1d + a2*cos(th1)*sin(th2)*th2d)^2 + a2^2*cos(th2)^2*th2d^2))/2 + (a1^2*m1*th1d^2)/2 - a2*g*m1*sin(th2);
R = 0; % Friction term
par ={ }; % System parameters
Q_i = {0 0}; Q_e = {0 0};
X={th1 th1d th2 th2d};
VF = EulerLagrange(L,X,Q_i,Q_e,R,par,'m','s','ode1');
pause(10) %Delay so the code is generated
tspan = [0 10];
y0=[10*pi/180;0;40*pi/180;0];
[t,Y] = ode113(@ode1, tspan, y0);
subplot(2,2,1);plot (t,Y(:,1)*180/pi);
subplot(2,2,2);plot (t,Y(:,2)*180/pi);
subplot(2,2,3);plot (t,Y(:,3)*180/pi);
subplot(2,2,4);plot (t,Y(:,4)*180/pi);
