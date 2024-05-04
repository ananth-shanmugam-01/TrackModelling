function [modelFunction, states, controls, parameters] = trackDynamics()

    Curv    = casadi.SX.sym('Curv');   % Curvature;
    theta   = casadi.SX.sym('theta');  % Racing line tangent angle
    xi      = casadi.SX.sym('xi');     % Track X-position - cartesian coordinates
    yi      = casadi.SX.sym('yi');     % Track Y-position - cartesian coordinates
    states.syms  = [Curv; theta; xi; yi];
    states.names = {'Curv'; 'theta'; 'xi'; 'yi'};
    states.lb    = [-3; -3*pi; -inf; -inf];
    states.ub    = [3; 3*pi; inf; inf];
    states.scale = [1; 0.1; 0.01; 0.01];
    states.num_x = numel(states.names);
    
    % Controls
    u                = casadi.SX.sym('u');      % Curvature Smoothing Factor
    controls.syms    = [u];
    controls.names   = {'u'};
    controls.lb    = [-inf];
    controls.ub    = [inf];
    controls.scale = [1];
    controls.num_u   = numel(controls.names);
    
    % Parameters
    xc          = casadi.SX.sym('xc');
    yc          = casadi.SX.sym('yc');
    c           = casadi.SX.sym('c'); % smoothing penalty
    parameters.syms  = [xc; yc ; c];
    parameters.names = {'xc', 'yc', 'c'};
    parameters.num_g = numel(parameters.names);

    % Algebraic Equations

    % Define system dynamics
    modelRhs = casadi.SX.sym('modelRhs',4);
    modelRhs(1) = u;
    modelRhs(2) = Curv;
    modelRhs(3) = cos(theta);
    modelRhs(4) = sin(theta);
    
    % Penalty Term
    modelPenalty = (xc - xi)^2 + (yc - yi)^2 + c*u^2;
    
    % Create Function to provide derivative and costs at each point
    modelFunction = casadi.Function('modelFunction',{states.syms,controls.syms,parameters.syms},{modelRhs, modelPenalty});


end