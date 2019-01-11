function [ xdot ] = pendulumStateDerivative(t, x , params )
% State derivative of pendulum model without short-range stiffness and
% reflexes
% x = [q qdot]'

I = params.I;
m = params.mass;
l = params.lc;
g = 9.81; % gravity

% Passive force to prevent joints from going outside their range of motion.
knee_r_range = params.knee_r_range;
klim = params.klim;
theta_ref = (knee_r_range(2) - knee_r_range(1))/2;
q = x(1) - (knee_r_range(1) + knee_r_range(2))/2;
Tlimit = -exp(klim*(q-theta_ref)) + exp(klim*(-q-theta_ref));

% Linear damping
Td = - params.d * x(2);

% Pendulum dynamics
xdot = zeros(2,1);
xdot(1) = x(2);
xdot(2) = 1/I *(-m*g*l*cos(x(1)) + Tlimit + Td + params.Tb);

end

