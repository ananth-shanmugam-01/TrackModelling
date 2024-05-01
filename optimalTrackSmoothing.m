%% rocket landing - Volkl, Muehlmeier
clear; clc;
import casadi.*

%% Load Track Data File
load('C:\Users\admin\Desktop\GitHub\QSS-LTS\dataFiles\230722_Endurance_lap2.mat')

lat_pos     = smoothdata(data_lap2.GPS_Latitude_deg,"gaussian","SmoothingFactor",0.001);
long_pos    = smoothdata(data_lap2.GPS_Longitude_deg,"gaussian","SmoothingFactor",0.001);
alt_pos     = data_lap2.GPS_Altitude_m;   

lat_origin  = lat_pos(1);
long_origin = long_pos(1);
alt_origin  = alt_pos(1);

origin = [lat_origin, long_origin, alt_origin];

measData = struct;

[measData.x, measData.y, measData.z] = latlon2local(lat_pos, long_pos, alt_pos, origin);

% sLap calculation from absolute distance
measData.dS = sqrt(gradient(measData.x).^2 + gradient(measData.y).^2);
measData.sLap = cumtrapz(measData.dS);

% Curvature calculation
x_dot   = gradient(measData.x)./gradient(measData.sLap);
x_ddot  = gradient(x_dot);
y_dot   = gradient(measData.y)./gradient(measData.sLap);
y_ddot  = gradient(y_dot);

measData.Curv = (x_dot.* y_ddot - x_ddot.*y_dot)./( x_dot.^2 + y_dot.^2).^1.5; 

% Trapezoidal Integration for track heading angle

measData.theta = cumtrapz(measData.sLap, measData.Curv);
measData.heading_origin = atan(y_dot(1)/x_dot(1));

measData.theta = measData.theta - measData.heading_origin;

%% Model Dynamics

% Declare model decision variables
% States
Curv    = SX.sym('Curv');   % Curvature;
theta   = SX.sym('theta');  % Racing line tangent angle
xi      = SX.sym('xi');     % Track X-position - cartesian coordinates
yi      = SX.sym('yi');     % Track Y-position - cartesian coordinates
states  = [Curv; theta; xi; yi];
stateNames = {'Curv'; 'theta'; 'xi'; 'yi'};
num_x = numel(states);

% Controls
u           = SX.sym('u');      % Curvature Smoothing Factor
controls    = [u];
controlNames = {'u'};
num_u = numel(controls);

% Parameters
xc           = SX.sym('xc');
yc           = SX.sym('yc');
c             = SX.sym('c'); % smoothing penalty
parameters  = [xc; yc ; c];
num_g       = numel(parameters);

% Define system dynamics
rhs = casadi.SX.sym('rhs',4);
rhs(1) = u;
rhs(2) = Curv;
rhs(3) = cos(theta);
rhs(4) = sin(theta);

% Penalty Term
cost = (xc - xi)^2 + (yc - yi)^2 + c*u^2;

% Create Function to provide derivative and costs at each point
f = casadi.Function('f',{states,controls,parameters},{rhs, cost});

%% Using Legendre Lagrangian Polynomials

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau = collocation_points(d, 'legendre');

% Collocation linear maps
[C,D,B] = collocation_coeff(tau);

%% Mesh Discretisation

% distance horizon, Control discretization
s = measData.sLap(end);        % m
N = 300;              % number of elements - phases are within this element
h = s/N;              % mesh size

% Initial Element Discretisation
% Constant spacing used for this case
dist = 0:h:s;

remesh = zeros(N,d+1);

for i = 1:numel(dist)-1
    dMesh = dist(i+1) - dist(i);
    remesh(i, 1) = dist(i);
    for i_ = 1:numel(tau)
        remesh(i, i_+1) = remesh(i, 1) + dMesh*tau(i_);
    end
end

dist_mesh = [reshape(remesh',1,[]), dist(end)];


%% NLP Formulation

% Start with an empty NLP
opti = Opti();

J = 0; % Objective term

% "Lift" initial conditions
Xk = opti.variable(num_x);
Uk = opti.variable(num_u);
Gk = opti.parameter(num_g);

% Collect all states/controls
Xs = {Xk};
Us = {Uk};
Gs = {Gk};

% Formulate the NLP
for k=1:N
  
   % Decision variables for helper states at each collocation point
   Xc = opti.variable(num_x, d);
   Uc = opti.variable(num_u, d);

   % Symbolic for parameter - not an decision variable
   Gc = opti.parameter(num_g,d);

   % Append these symbolics (decision variables + parameters)
   Xs{end+1} = Xc; % Controls
   Us{end+1} = Uc;
   Gs{end+1} = Gc;

   % Evaluate ODE right-hand-side at all helper states
   [ode, quad] = f(Xc, Uc, Gc);

   % No Additional Constraints Required 
   % opti.subject_to(f_Power(Xc, Uc) <= 1); 

   % Add contribution to quadrature function
   J = J + quad*B*h;

   % Get interpolating points of collocation polynomial
   Z_s = [Xk Xc]; % States
   Z_u = [Uk Uc]; % Controls

   % Get slope of interpolating polynomial (normalized)
   Pidot = (1/h)*Z_s*C;
   opti.subject_to(Pidot == ode); % Match with ODE right-hand-side 

   % State and control interpolated to end of collocation interval
   Xk_end = Z_s*D;
   Uk_end = Z_u*D;

   % New decision variable for state and control at end of interval
   Xk = opti.variable(num_x);
   Uk = opti.variable(num_u);
   Gk = opti.parameter(num_g);  
   
   Xs{end+1} = Xk;
   Us{end+1} = Uk;
   Gs{end+1} = Gk;

   opti.subject_to(Xk_end==Xk); % Continuity constraints
   opti.subject_to(Uk_end==Uk);

end

Xs = [Xs{:}];
Us = [Us{:}];
Gs = [Gs{:}];

%% Set Parameters
track_x     = interp1(measData.sLap, measData.x, dist_mesh,'pchip');
track_y     = interp1(measData.sLap, measData.y, dist_mesh,'pchip');
track_curv  = interp1(measData.sLap, measData.Curv, dist_mesh,'pchip');
track_theta = interp1(measData.sLap, measData.theta, dist_mesh, 'pchip');

smooth_weight   = 0.002;

opti.set_value(Gs(1,:), track_x);
opti.set_value(Gs(2,:), track_y);
opti.set_value(Gs(3,:), smooth_weight);

% Path Constraints
opti.subject_to(-0.5 <= Xs(1,:) <= 0.5);
opti.subject_to(-3*pi <= Xs(2,:) <= 3*pi);

% Initial & Terminal Constraints
% States
stateNames = {'Curv'; 'theta'; 'xi'; 'yi'};
opti.subject_to(Xs(1,1)     == Xs(1,end)); 
opti.subject_to(Xs(2,1)     == Xs(2,end)); 
opti.subject_to(Xs(3,1)     == track_x(1)); 
opti.subject_to(Xs(3,end)   == track_x(end)); 
opti.subject_to(Xs(4,1)     == track_y(1));
opti.subject_to(Xs(4,end)   == track_y(end));

% Controls

% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial Conditions
opti.set_initial(Xs(1,:), track_curv  .*3);
opti.set_initial(Xs(2,:), track_theta .*2);
opti.set_initial(Xs(3,:), track_x);
opti.set_initial(Xs(4,:), track_y);


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Objective
opti.minimize(J);
plugin_opts = struct();
solver_opts = struct('max_iter',500);
opti.solver('ipopt',plugin_opts, solver_opts);

sol = opti.solve();

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Create Results Struct

results = struct;

for i = 1:num_x
    results.states.(stateNames{i}) = opti.value(Xs(i,:));
end

for i = 1:num_u
    results.controls.(controlNames{i}) = opti.value(Us(i,:));
end

results.sLap   = dist_mesh;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot results

figure(1); clf
hold on
scatter(track_x, track_y,'DisplayName','Input')
scatter(results.states.xi, results.states.yi,'DisplayName','OC')
hold off; legend
xlabel('[m]')
ylabel('[m]')
title('Track Map')
grid minor; box on;

figure(2); clf
hold on
plot(measData.sLap, measData.Curv,'DisplayName','Input')
plot(dist_mesh, track_curv,'DisplayName','Initial Guess')
plot(results.sLap, results.states.Curv,'DisplayName','OC')
hold off; legend
title('Track Curvature')
ylabel('[1/m]')
grid minor; box on;

figure(3); clf
hold on
plot(measData.sLap, measData.theta,'DisplayName','Input')
plot(dist_mesh, track_theta,'DisplayName','Initial Guess')
plot(results.sLap, results.states.theta,'DisplayName','OC')
hold off; legend
title('Track Heading Angle')
ylabel('[rad]')
grid minor; box on;