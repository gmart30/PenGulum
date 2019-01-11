function [ xdot ] = pendulumStateDerivative_SRS_Ffb(t, x , Z, dZ)
% State derivative of pendulum model with short-range stiffness and
% force-related reflexes
% x = [q qdot]'
% x = [q qdot TSRS Tref]

global hist 
global params

I = params.I;
m = params.mass;
l = params.lc;
g = 9.81; % gravity

xdot = zeros(4,1);

% Passive force to prevent joints from goint outside their range of motion.
theta_ref = (params.knee_r_range(2) - params.knee_r_range(1))/2;
q = x(1) - (params.knee_r_range(1) + params.knee_r_range(2))/2;
Tlimit = -exp(params.klim*(q-theta_ref)) + exp(params.klim*(-q-theta_ref));

% Linear damping
Td = - params.d * x(2);
Td_delay = - params.d * Z(2); % Td(t - tau) with tau reflex delay
Td_der_delay = - params.d * dZ(2); % Td_dot(t - tau) with tau reflex delay

% SRS
kSRS = params.kSRS * params.Tb; 
delta_theta = x(1) - params.theta0;

s_velneg = 1/2*(erf(-10*x(2)) + 1);
if hist == 1 & s_velneg > 0.001 % Short-range stiffness only during first swing excursion
    s_stcrit = 1/2*(erf(50*(delta_theta + params.delta_theta_crit)) + 1);
    s_gtcrit = 1/2*(erf(50*(-delta_theta- params.delta_theta_crit)) + 1);    
    xdot_SRS = 1/0.01 * ( s_stcrit*(-x(3) - kSRS * delta_theta) + s_gtcrit*(-x(3) + kSRS * params.delta_theta_crit));
    
    s_velpos = 1/2*(erf(10*x(2)) + 1);
    xdot(3) = s_velneg * xdot_SRS + s_velpos * 1/params.tauSRS * (-x(3));   
else
    xdot(3) = 1/params.tauSRS * (-x(3));
    hist = 0;
end

u = max(0, params.kF * (Z(3) + Td_delay) + params.kdF * (dZ(3)+Td_der_delay)); % torque input

xdot(1) = x(2);
xdot(2) = 1/I *(-m*g*l*cos(x(1)) + Tlimit + x(3) + x(4) + Td + params.Tb); % pendulum dynamics
xdot(4) = (u - x(4))/0.01; % derivative of reflex torque

end

