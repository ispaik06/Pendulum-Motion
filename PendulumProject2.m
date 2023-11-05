%% Double Pendulum
clear; close all; clc;

% parameter settings
gVal = 9.81;
m1 = 1;
m2 = 1;

l1 = 1;
l2 = 1;

theta1 = pi*95/180;
theta2 = pi/4;

theta1Vel = 0;
theta2Vel = 0;

count = 1400;

time_step = 20/count;
time_span = linspace(0, 20, count);

theta1_List = zeros(1, count+1);
theta2_List = zeros(1, count+1);
theta1_List(1) = theta1;
theta2_List(1) = theta2;

% Numerical Calculation
idx = 2;
for i=linspace(0, 10, count)
    [theta1, theta1Vel] = theta1_cal(m1, m2, l1, l2, gVal, theta1, theta2, theta1Vel, theta2Vel, time_step);
    [theta2, theta2Vel] = theta2_cal(m1, m2, l1, l2, gVal, theta1, theta2, theta1Vel, theta2Vel, time_step);

    theta1_List(idx) = theta1;
    theta2_List(idx) = theta2;
    idx = idx+1;
end

x1 = l1.*sin(theta1_List);
y1 = -l1.*cos(theta1_List);
x2 = x1 + l2.*sin(theta2_List);
y2 = y1 - l2.*cos(theta2_List);

%% Motion Visualization
close all;

figure;

for i=1:count+1
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r');
    hold on;
    plot([0, x1(i)], [0, y1(i)], 'r-', 'LineWidth', 2);
    plot([x1(i), x2(i)], [y1(i), y2(i)], 'r-', 'LineWidth', 2);
    plot(x1(i), y1(i), 'yo', 'MarkerFaceColor', 'y');
    plot(x2(i), y2(i), 'go', 'MarkerFaceColor', 'g');
    set(gca, 'Color', 'k');
    title('\bfDouble Pendulum Motion');
    axis([-3 3 -3 3]);
    axis square;
    pause(0.001);
    clf;
end

%% Trajectory
figure;
for i=1:count+1
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r');
    hold on;
    plot(x1(i), y1(i), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 1);
    plot(x2(i), y2(i), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 1);
    set(gca, 'Color', 'k');
    title('\bfDouble Pendulum Trajectory');
    axis([-3 3 -3 3]);
    axis square;
    pause(0.001);
end



%% FUNCTIONS CODE
function theta1_ddot = theta1_acceleration(m1, m2, l1, l2, g, theta1, theta2, theta_t1, theta_t2)
theta1_ddot = -(g*m1*sin(theta1) + g*m2*sin(theta1) - g*m2*cos(theta1 - theta2)*sin(theta2) + l2*m2*theta_t2^2*sin(theta1 - theta2) + l1*m2*theta_t1^2*cos(theta1 - theta2)*sin(theta1 - theta2))/(l1*(m1 + m2 - m2*cos(theta1 - theta2)^2));
end


function theta2_ddot = theta2_acceleration(m1, m2, l1, l2, g, theta1, theta2, theta_t1, theta_t2)
theta2_ddot = (g*m1*cos(theta1 - theta2)*sin(theta1) - g*m2*sin(theta2) - g*m1*sin(theta2) + g*m2*cos(theta1 - theta2)*sin(theta1) + l1*m1*theta_t1^2*sin(theta1 - theta2) + l1*m2*theta_t1^2*sin(theta1 - theta2) + l2*m2*theta_t2^2*cos(theta1 - theta2)*sin(theta1 - theta2))/(l2*(m1 + m2 - m2*cos(theta1 - theta2)^2));
end


function [theta1, theta_t1] = theta1_cal(m1, m2, l1, l2, g, theta1, theta2, theta_t1, theta_t2, timeStep)
theta_t1 = theta_t1 + timeStep*theta1_acceleration(m1,m2,l1,l2,g,theta1,theta2,theta_t1,theta_t2);
theta1 = theta1 + timeStep*theta_t1;
end

function [theta2, theta_t2] = theta2_cal(m1, m2, l1, l2, g, theta1, theta2, theta_t1, theta_t2, timeStep)
theta_t2 = theta_t2 + timeStep*theta2_acceleration(m1,m2,l1,l2,g,theta1,theta2,theta_t1,theta_t2);
theta2 = theta2 + timeStep*theta_t2;
end
