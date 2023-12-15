%% Quiz4 Solution
clc
clear all
close all

syms th1(t) th2(t) th1d th2d t a1 a2 m1 m2 g

%Part A
% Rotation Matrices
R01=[cos(th1) -sin(th1) 0;sin(th1) cos(th1) 0;0 0 1]*[1 0 0;0 0 -1;0 1 0];
R12=[cos(th2) -sin(th2) 0;sin(th2) cos(th2) 0;0 0 1]*[1 0 0;0 1 0;0 0 1];

%Displacement Vectors
D01=[a1*cos(th1);a1*sin(th1);0];
D12=[a2*cos(th2);a2*sin(th2);0];

R02=R01*R12;

%Homogenous Transformation Matrices 
add=[0 0 0 1];

H01=[R01 D01;add];
H12=[R12 D12;add];

H02temp=H01*H12;
H02=H02temp(t);
D02=H02(1:3,4)% Displacement Vector of End Effector

%% Part B

disp_dot = diff(D02,t);% Differentiation of the displacement vector
disp_dot=subs(disp_dot,[diff(th1(t),t), diff(th2(t),t)],[sym('th1d'), sym('th2d')]);

linjacobian=jacobian(disp_dot,[th1d,th2d])% Jacobian Matrix of Linear Velocity

R00=eye(3);
omega_vec=R00*[0;0;th1d]+R01*[0;0;th2d];
angjacobian=jacobian(omega_vec,[th1d,th2d])% Jacobian Matrix of Angular Velocity

%% Part C
th1const = linspace(0,pi/2,100);
th2const = linspace(pi/2,-pi/2,100);
a1const = 5;
a2const = 3;
D02const = [a1const.*cos(th1const)+a2const.*cos(th1const).*cos(th2const);a1const.*sin(th1const)+a2const.*cos(th2const).*sin(th1const);a2const.*sin(th2const)];
figure(1)
plot3(D02const(1,:),D02const(2,:),D02const(3,:),'r-', 'LineWidth', 4)
hold on
grid on
title('End Effector Trajectory');
xlabel('X02 variation');
ylabel('Y02 variation');
zlabel('Z02 variation');
hold off

%% Part E
syms th1d(t) th2d(t) th1dd th2dd
K1=0.5*m1*(a1^2)*(th1d(t).^2);
U1=0;

V2_2=simplify(disp_dot(1,1).^2+disp_dot(2,1).^2+disp_dot(3,1).^2);% used the ouptput of the expression in next line
V2_2=simplify((a1*th1d(t)*sin(th1(t)) + a2*th1d(t)*cos(th2(t))*sin(th1(t)) + a2*th2d(t)*cos(th1(t))*sin(th2(t)))^2 + (a1*th1d(t)*cos(th1(t)) + a2*th1d(t)*cos(th1(t))*cos(th2(t)) - a2*th2d(t)*sin(th1(t))*sin(th2(t)))^2 + a2^2*th2d(t)^2*cos(th2(t))^2);
K2=0.5*m2*V2_2;
U2=m1*g*D02(3,1);

L=simplify((K1+K2)-(U1+U2))% Lagrangian expression

%Partial derivatives wrt 'derivative of joint variable'
dL_dth1d=functionalDerivative(L,th1d);
dL_dth2d=functionalDerivative(L,th2d);

%Partial derivatives wrt 'joint variable'
dL_dth1=functionalDerivative(L,th1);
dL_dth2=functionalDerivative(L,th2);

%Derivatives with respect to time
ddL_dth1ddt=diff(dL_dth1d,t);
ddL_dth1ddt=simplify(subs(ddL_dth1ddt,[diff(th1(t),t), diff(th2(t),t), diff(th1d(t),t)],[sym('th1d'), sym('th2d'), sym('th1dd')]));

ddL_dth2ddt=diff(dL_dth2d,t);
ddL_dth2ddt=simplify(subs(ddL_dth2ddt,[diff(th1(t),t), diff(th2(t),t), diff(th2d(t),t)],[sym('th1d'), sym('th2d'), sym('th2dd')]));

% Equations of motion
Tau_1=simplify(ddL_dth1ddt-dL_dth1) %Motion for Revolute Joint 1

Tau_2=simplify(ddL_dth2ddt-dL_dth2) %Motion for Revolute Joint 2








