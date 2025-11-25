%% Optimal Smoothing Track Model

clear; clc;
import casadi.*

%% Load Track Data File

load("C:\Users\Ananth\Desktop\GitHub\DataFiles\Track\FS\230722_Endurance_lap2.mat")

trackCoordinates = struct;
trackCoordinates.latitude_deg   = data_lap2.GPS_Latitude_deg;
trackCoordinates.longitude_deg  = data_lap2.GPS_Longitude_deg;
trackCoordinates.altitude_m     = data_lap2.GPS_Altitude_m;

measData = initialiseTrackData(trackCoordinates);

%% Model Function, Penalty, States, Controls, Parameters

[modelFunction, states, controls, parameters] = trackDynamics();

%% Create Lagrange Polynomials

collocation = struct;

% Degree of interpolating polynomial
collocation.degree = 3;

% Get collocation points
collocation.tau = collocation_points(collocation.degree, 'legendre');

% Collocation linear maps
[collocation.C,collocation.D,collocation.B] = collocation_coeff(collocation.tau);
% C - derivative coefficients of polynomial
% D - interpolation coefficients
% B - integral coefficients

%% Mesh Discretisation

mesh = struct;

% distance horizon, Control discretization
mesh.endPoint       = measData.sLap(end);               % m
mesh.numIntervals   = 500;                              % number of elements - phases are within this element
mesh.dMesh          = mesh.endPoint/mesh.numIntervals;  % Interval Size

% Initial Element Discretisation
% Constant spacing used for this case
mesh.meshIntervals = 0:mesh.dMesh:mesh.endPoint;

remesh = zeros(mesh.numIntervals,collocation.degree+1);

for i = 1:numel(mesh.meshIntervals)-1
    dMeshLocal = mesh.meshIntervals(i+1) - mesh.meshIntervals(i);
    remesh(i, 1) = mesh.meshIntervals(i);
    for i_ = 1:numel(collocation.tau)
        remesh(i, i_+1) = remesh(i, 1) + dMeshLocal*collocation.tau(i_);
    end
end

mesh.meshPoints = [reshape(remesh',1,[]), mesh.meshIntervals(end)];

clear dMeshLocal remesh i i_

%% Combine into Settings

settings = struct;
settings.collocation = collocation;
settings.mesh        = mesh;
settings.states      = states;
settings.controls    = controls;
settings.parameters  = parameters;

%% Penalty Weight 

smooth_weight   = 2;

%% Create Optimisation Problem
[opti, cost, Xs, Us, Gs] = createOptiProblem(settings, modelFunction, []);

%% Set Parameters - Constants

constants   = struct;
constants.track_x     = interp1(measData.sLap, measData.x,     mesh.meshPoints,'pchip');
constants.track_y     = interp1(measData.sLap, measData.y,     mesh.meshPoints,'pchip');
constants.track_curv  = interp1(measData.sLap, measData.Curv,  mesh.meshPoints,'pchip');
constants.track_theta = interp1(measData.sLap, measData.theta, mesh.meshPoints,'pchip');

opti.set_value(Gs(1,:), constants.track_x);
opti.set_value(Gs(2,:), constants.track_y);
opti.set_value(Gs(3,:), smooth_weight);

%% Initial & Terminal Constraints

opti.subject_to(Xs(1,1)     == Xs(1,end)); % Initial and Final Curvatures aligned
% opti.subject_to(Xs(2,1)     == Xs(2,end)); % Initial and Final Heading Angles aligned
opti.subject_to(Xs(3,1)     == constants.track_x(1)); 
opti.subject_to(Xs(3,end)   == constants.track_x(end)); 
opti.subject_to(Xs(4,1)     == constants.track_y(1));
opti.subject_to(Xs(4,end)   == constants.track_y(end));

%% Set Initial Solution

% Copied over to maintain idealised structure, where initial solution is
% provided externally and will be interpolated to be coherent with the OC
% mesh characteristics

initSolution = struct;
initSolution.track_x     = interp1(measData.sLap, measData.x,     mesh.meshPoints,'pchip');
initSolution.track_y     = interp1(measData.sLap, measData.y,     mesh.meshPoints,'pchip');
initSolution.track_curv  = interp1(measData.sLap, measData.Curv,  mesh.meshPoints,'pchip');
initSolution.track_theta = interp1(measData.sLap, measData.theta, mesh.meshPoints,'pchip');
 

opti.set_initial(Xs(1,:), initSolution.track_curv);
opti.set_initial(Xs(2,:), initSolution.track_theta);
opti.set_initial(Xs(3,:), initSolution.track_x);
opti.set_initial(Xs(4,:), initSolution.track_y);

%% IPOPT Solver Options
pluginOptions = struct;

solverOptions = struct;
solverOptions.max_iter = 500;

opti.solver('ipopt',pluginOptions, solverOptions);

%% Run Optimisation Problem
sol = opti.solve();

%% Create Results Struct

results = struct;
Track = struct; % Save Track Model Specific Outputs

for i = 1:states.num_x
    results.states.(states.names{i}) = opti.value(Xs(i,:));
    Track.(states.names{i})          = opti.value(Xs(i,:));
end

for i = 1:controls.num_u
    results.controls.(controls.names{i}) = opti.value(Us(i,:));
end

for i = 1:parameters.num_g
    results.parameters.(parameters.names{i}) = opti.value(Gs(i,:));
end

Track.sLap      = settings.mesh.meshPoints;

results.mesh = settings.mesh;
results.settings = settings;
results.initSolution = initSolution;

save('FSUK_2023.mat','Track')
%% Plot results

figure(1); clf
hold on
scatter(measData.x, measData.y,'DisplayName','Input')
scatter(results.states.xi, results.states.yi,'DisplayName','OC')
hold off; legend
xlabel('[m]')
ylabel('[m]')
title('Track Map')
grid minor; box on;

figure(2); clf
hold on
plot(measData.sLap, measData.Curv,'DisplayName','Input')
plot(mesh.meshPoints, initSolution.track_curv,'DisplayName','Initial Guess')
plot(results.mesh.meshPoints, results.states.Curv,'DisplayName','OC')
hold off; legend
title('Track Curvature')
ylabel('[1/m]')
grid minor; box on;

figure(3); clf
hold on
plot(measData.sLap, measData.theta,'DisplayName','Input')
plot(mesh.meshPoints, initSolution.track_theta,'DisplayName','Initial Guess')
plot(results.mesh.meshPoints, results.states.theta,'DisplayName','OC')
hold off; legend
title('Track Heading Angle')
ylabel('[rad]')
grid minor; box on;