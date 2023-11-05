%% Linear Pendulum
% 단진자 운동방정식에서 작은각 근사를 사용하여
% 선형 미분방정식을 구하고, dsolve로 해를 직접 구하는 방법

clear; close all; clc;

syms g l theta(t)
angVel= diff(theta);
eqn = diff(theta, 2) + g/l*sin(theta) == 0;

% 작은 값 근사
syms x
approx = taylor(sin(x), x, 'Order', 2);
approx = subs(approx, x, theta(t));
eqnAprroxed = subs(eqn, sin(theta), approx);

% 미분방정식 해 구하기
theta_0Val = pi/3; angVel_0 = 0;
cond = [theta(0)==theta_0Val, angVel(0)==angVel_0];
thetaSolved(t) = dsolve(eqnAprroxed, cond);

syms x(t) y(t)
x(t) = l*sin(thetaSolved(t));
y(t) = -l*cos(thetaSolved(t));

% 값 대입
gVal = 9.81; lVal = 1;
thetaReal(t) = subs(thetaSolved(t), [g, l], [gVal, lVal]);
xPos(t) = subs(x(t), [g, l], [gVal, lVal]);
yPos(t) = subs(y(t), [g, l], [gVal, lVal]);

omega_0Val = sqrt(gVal/lVal);
T0 = 2*pi/(omega_0Val);


figure;
ezplot(thetaReal(t)/pi*180, [0 5*T0]);
title('\theta - \itt \rm\bfgraph for linear pendulum')
xlabel('\itt(s)');
ylabel('\theta (°)');


% time = 0:0.1:5*T;
% thetaReal = subs(thetaReal(t), t, time);
% thetaDeg = double(thetaReal)/pi*180;
% 
% figure;
% plot(time, thetaDeg, '-r');
% title('\theta - \itt \rm\bfgraph for linear pendulum')
% xlabel('\itt (s)');
% ylabel('\theta (°)');

%% Motion Visualization
close all;
figure;
for i=0:0.1:5*T0
    plot([-0.7 0.7], [0 0], 'b-', 'LineWidth', 5);
    hold on;
    plot([0 xPos(i)], [0 yPos(i)], 'r-', 'LineWidth', 2);
    plot(xPos(i), yPos(i), 'go', 'MarkerFaceColor', 'g', 'Markersize', 13);
    set(gca, 'color', 'k')
    title('\bfLinear Pendulum Motion')
    xlabel('\itx\rm(m)')
    ylabel('\ity\rm(m)')
    axis([-2 2 -2 1]);
    axis square;
    pause(0.001);
    clf;
end

close all;


%% NonLinear Pendulum 1
% 단진자 운동방정식(비선형 미분방정식)의 이론적인 해 식을 이용
% 'Exact solution for the nonlinear pendulum' 논문 참고

close all; clc;

syms theta_0 omega_0
syms theta_Theory(t)

theta_Theory(t) = 2*asin(sin(theta_0/2)* ellipj(ellipke(sin(theta_0/2)^2)-omega_0*t, sin(theta_0/2)^2) );

syms x_N(t) y_N(t)
x_N(t) = l*sin(theta_Theory(t));
y_N(t) = -l*cos(theta_Theory(t));


thetaTheory(t) = subs(theta_Theory(t), [g, l, theta_0, omega_0], [gVal, lVal, theta_0Val, omega_0Val]);
xPos_N(t) = subs(x(t), [g, l, theta_0, omega_0], [gVal, lVal, theta_0Val, omega_0Val]);
yPos_N(t) = subs(y(t), [g, l, theta_0, omega_0], [gVal, lVal, theta_0Val, omega_0Val]);

T = (2/pi)*T0*ellipke(sin(theta_0Val/2)^2);

figure;
ezplot(thetaReal(t)/pi*180, [0 5*T0]);
hold on;
ezplot(thetaTheory(t)/pi*180, [0 5*T]); % theta_0Val 값이 작을수록 그래프가 겹침
title('\theta - \itt \rm\bfgraph')
xlabel('\itt(s)');
ylabel('\theta (°)');
legend('Linear', 'NonLinear');

%% Motion Visualization
for i=0:0.1:5*T
    plot([-0.7 0.7], [0 0], 'b-', 'LineWidth', 5);
    hold on;
    plot([0 xPos(i)], [0 yPos(i)], 'r-', 'LineWidth',2);
    plot(xPos(i), yPos(i), 'go', 'MarkerFaceColor', 'g', 'Markersize', 13);
    set(gca, 'color', 'k')
    title('\bfNonlinear Pendulum Motion')
    xlabel('\itx\rm(m)')
    ylabel('\ity\rm(m)')
    axis([-2 2 -2 1]);
    axis square;
    pause(0.001);
    clf;
end

close all;


%% NonLinear Pendulum 2
clear; close all; clc;

syms theta theta_dot omega_0
gVal = 9.81; lVal = 1;

omega_0Val = sqrt(gVal/lVal);
E(theta, theta_dot, omega_0) = (1/2)*(theta_dot^2+(2*omega_0*sin(theta/2))^2);
E(theta, theta_dot) = subs(E,omega_0,omega_0Val);

figure;
fcontour(E(pi*theta, omega_0Val*theta_dot), 2*[-2 2 -2 2], 'LineWidth', 1);
grid on;
title('Constant Energy Contours in Phase Space');
xlabel('\theta/\pi');
ylabel('\theta_{dot}/\omega_0');


%% 구현
% 2차 미분방정식인 비선형 운동방정식을 ode45 solver를 이용해 수치적으로 푼다.
clear; close all; clc;

syms m g l 
syms theta(t) theta_t(t) omega_0
eq = [theta_t == diff(theta),  diff(theta_t) == -omega_0^2*sin(theta)];

% 값 대입
gVal = 9.81; lVal = 1; 
omega_0Val = sqrt(gVal/lVal);

eq = subs(eq, omega_0, omega_0Val);
vars = [theta, theta_t];

[M, F] = massMatrixForm(eq,vars);
f = M\F;
f = odeFunction(f, vars);

%% 1. Closed Energy Contours에 대한 운동방정식 풀기
cond = [0; 1.993*omega_0Val];
t_start = 0;
t_end = 10;
solved = ode45(f, [t_start, t_end], cond);

% Closed Path in Phase Space 그리기
figure(1);
yyaxis left;
plot(solved.x, solved.y(1,:), '-ro');
ylabel('\theta (rad)', 'Color','r');

yyaxis right;
plot(solved.x, solved.y(2,:), '-ko');
ylabel('\theta_t (rad/s)', 'Color', 'k');

title('\bfClosed Path in Phase Space');
xlabel('\itt (s)');
grid on;

%% Motion Visualization (Closed Path)
xPos = @(t) lVal*sin(deval(solved, t, 1));
yPos = @(t) -lVal*cos(deval(solved, t, 1));

close all;
figure(2);
for i=t_start:0.1:t_end
    plot([-0.7 0.7], [0 0], 'b-', 'LineWidth', 5);
    hold on;
    plot([0 xPos(i)], [0 yPos(i)], 'r-', 'LineWidth',2);
    plot(xPos(i), yPos(i), 'go', 'MarkerFaceColor', 'g', 'Markersize', 13);
    set(gca, 'color', 'k')
    title('\bfNonlinear Pendulum Motion (Closed Path)')
    xlabel('\itx\rm(m)')
    ylabel('\ity\rm(m)')
    axis([-2 2 -2 2]);
    axis square;
    pause(0.001);
    clf;
end

close all;

%% 2. Open Energy Contours에 대한 운동방정식 풀기
cond = [0; 2.001*omega_0Val];
t_start = 0;
t_end = 10;
solved = ode45(f, [t_start, t_end], cond);

% Open Path in Phase Space 그리기
figure(1);
yyaxis left;
plot(solved.x, solved.y(1,:), '-ro');
ylabel('\theta (rad)', 'Color','r');

yyaxis right;
plot(solved.x, solved.y(2,:), '-ko');
ylabel('\theta_t (rad/s)', 'Color','k');

grid on;
title('\bfOpen Path in Phase Space');
xlabel('\itt (s)');


%% Motion Visualization (Open Path)
xPos = @(t) lVal*sin(deval(solved, t, 1));
yPos = @(t) -lVal*cos(deval(solved, t, 1));

close all;
figure(2);
for i=t_start:0.1:t_end
    plot([-0.7 0.7], [0 0], 'b-', 'LineWidth', 5);
    hold on;
    plot([0 xPos(i)], [0 yPos(i)], 'r-', 'LineWidth',2);
    plot(xPos(i), yPos(i), 'go', 'MarkerFaceColor', 'g', 'Markersize', 13);
    set(gca, 'color', 'k')
    title('\bfNonlinear Pendulum Motion (Open Path)')
    xlabel('\itx\rm(m)')
    ylabel('\ity\rm(m)')
    axis([-2 2 -2 2]);
    axis square;
    pause(0.001);
    clf;
end

close all;
