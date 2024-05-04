function [opti, cost, Xs, Us, Gs] = createOptiProblem(settings, modelFunction, constraintFunction)

    % NLP Formulation
    
    % Start with an empty NLP
    opti = casadi.Opti();
    
    cost = 0; % Objective term
    
    % "Lift" initial conditions
    Xk = opti.variable(settings.states.num_x);
    Uk = opti.variable(settings.controls.num_u);
    Gk = opti.parameter(settings.parameters.num_g);
    
    % Collect all states/controls
    Xs = {Xk};
    Us = {Uk};
    Gs = {Gk};
    
    % Formulate the NLP
    for k=1:settings.mesh.numIntervals  
      
       % Decision variables for helper states at each collocation point
       Xc = opti.variable(settings.states.num_x, settings.collocation.degree);
       Uc = opti.variable(settings.controls.num_u, settings.collocation.degree);
    
       % Symbolic for parameter - not an decision variable
       Gc = opti.parameter(settings.parameters.num_g, settings.collocation.degree);
    
       % Append these symbolics (decision variables + parameters)
       Xs{end+1} = Xc; % Controls
       Us{end+1} = Uc;
       Gs{end+1} = Gc;
    
       % Evaluate ODE right-hand-side at all helper states
       [ode, quad] = modelFunction(Xc, Uc, Gc);
    
% %        Improve this formulation to read additional functions
%        if exist("constraintFunction")
%            opti.subject_to(constraintFunction(Xc, Uc, Gc) <= 1);
%        end

       % Add contribution to quadrature function
       cost = cost + quad * settings.collocation.B * settings.mesh.dMesh ;
    
       % Get interpolating points of collocation polynomial
       Z_s = [Xk Xc]; % States
       Z_u = [Uk Uc]; % Controls
    
       % Get slope of interpolating polynomial (normalized)
       Pidot = (1/settings.mesh.dMesh ) * Z_s * settings.collocation.C;
       opti.subject_to(Pidot == ode); % Match with ODE right-hand-side 
    
       % State and control interpolated to end of collocation interval
       Xk_end = Z_s * settings.collocation.D;
       Uk_end = Z_u * settings.collocation.D;
    
       % New decision variable for state and control at end of interval
       Xk = opti.variable(settings.states.num_x);
       Uk = opti.variable(settings.controls.num_u);
       Gk = opti.parameter(settings.parameters.num_g);  
       
       Xs{end+1} = Xk;
       Us{end+1} = Uk;
       Gs{end+1} = Gk;
    
       opti.subject_to(Xk_end==Xk); % Continuity constraints
       opti.subject_to(Uk_end==Uk);
    
    end
    
    Xs = [Xs{:}];
    Us = [Us{:}];
    Gs = [Gs{:}];

    % Decision Variable Bounds
    for i = 1:numel(settings.states.num_x)
        opti.subject_to(settings.states.lb(i) <= Xs(i,:) <= settings.states.ub(i));
    end

    for i = 1:numel(settings.controls.num_u)
        opti.subject_to(settings.controls.lb(i) <= Us(i,:) <= settings.controls.ub(i));
    end
    
    % Objective
    opti.minimize(cost);

end
