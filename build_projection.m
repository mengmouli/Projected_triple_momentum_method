function proj_fun = build_projection(constraint_fun)

    if isempty(constraint_fun)
        % No constraints — identity projection
        proj_fun = @(x) x;
        return;
    end

    proj_fun = @(x) project_yalmip(x,constraint_fun);
end
