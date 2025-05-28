function x_proj = project_yalmip(x_val, constraint_fun)
    d = length(x_val);
    y = sdpvar(d, 1);  % create optimization variable

    Constraints = [];
    c_list = constraint_fun(y);
    Constraints = [Constraints, c_list{:}];
    
    Objective = (y - x_val)' * (y - x_val);
    % assign(y, x_val);
    
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);
    ops.mosek.MSK_DPAR_INTPNT_TOL_REL_GAP    = 1e-10;
    ops.mosek.MSK_DPAR_INTPNT_TOL_PFEAS      = 1e-10;
    ops.mosek.MSK_DPAR_INTPNT_TOL_DFEAS      = 1e-10;
    ops.mosek.MSK_DPAR_INTPNT_TOL_MU_RED     = 1e-10;
    ops.mosek.MSK_DPAR_INTPNT_TOL_PATH       = 1e-10;
    ops.mosek.MSK_DPAR_INTPNT_TOL_STEP_SIZE  = 1e-10;

    % ops = sdpsettings('solver', 'sedumi', 'verbose', 0, 'sedumi.eps', 1e-12);

    sol = optimize(Constraints, Objective, ops);
    if sol.problem ~= 0
        warning(['Projection failed: ', sol.info]);
    end

    x_proj = value(y);
end
