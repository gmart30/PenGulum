function [ xdot ] = pendulumStateDerivative_SRS(t, x , params )
% State derivative of pendulum model with short-range stiffness but without
% reflexes
% x = [q qdot]'
% x = [q qdot TSRS]'
global hist

I = params.I;
m = params.mass;
l = params.lc;
g = 9.81; % gravity

xdot = zeros(3,1);

% Passive force to prevent joints from going outside their range of motion.
knee_r_range = params.knee_r_range;
klim = params.klim;
theta_ref = (knee_r_range(2) - knee_r_range(1))/2;
q = x(1) - (knee_r_range(1) + knee_r_range(2))/2;
Tlimit = -exp(klim*(q-theta_ref)) + exp(klim*(-q-theta_ref));

% Linear damping
Td = - params.d * x(2);

% Pendulum dynamics
xdot(1) = x(2);
xdot(2) = 1/I *(-m*g*l*cos(x(1)) + Tlimit + x(3) + Td + params.Tb);

% SRS
kSRS = params.kSRS * params.Tb; 
delta_theta = x(1) - params.theta0;

s_velneg = 1/2*(erf(-10*x(2)) + 1);
if hist == 1 & s_velneg > 0.001 % Short-range stiffness only during first swing excursion
    s_stcrit = 1/2*(erf(50*(delta_theta + params.delta_theta_crit)) + 1);
    s_gtcrit = 1/2*(erf(50*(-delta_theta- params.delta_theta_crit)) + 1);    
    xdot_SRS = 1/0.001 * ( s_stcrit*(-x(3) - kSRS * delta_theta) + s_gtcrit*(-x(3) + kSRS * params.delta_theta_crit));
    s_velpos = 1/2*(erf(10*x(2)) + 1);
    xdot(3) = s_velneg * xdot_SRS + s_velpos * 1/params.tauSRS * (-x(3));   
else
    xdot(3) = 1/params.tauSRS * (-x(3));
    hist = 0;
end

end

