syms m1 m2 l1 l2 g theta1 theta2 theta_t1 theta_t2 theta_tt1 theta_tt2

e1 = (m1+m2)*l1^2*theta_tt1 + m2*l1*l2*theta_tt2*cos(theta1-theta2) + m2*l1*l2*theta_t2^2*sin(theta1-theta2) + (m1+m2)*g*l1*sin(theta1) == 0;
e2 = m2*l2^2*theta_tt2 + m2*l1*l2*theta_tt1*cos(theta1-theta2) - m2*l1*l2*theta_t1^2*sin(theta1-theta2) + m2*g*l2*sin(theta2) ==0 ;

sol = solve([e1, e2], [theta_tt1, theta_tt2]);

theta_tt1_sol = solutions.theta_tt1;
theta_tt2_sol = solutions.theta_tt2;


%% Coupled Pendulum


%% Sperical Pendulum