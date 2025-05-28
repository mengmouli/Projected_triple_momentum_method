function results = proj_gradient(A, B, C, D, x0, num_steps, grad_f, proj_fun)
    d = length(x0);
    % Kronecker
    Ad = kron(A, eye(d));
    Bd = kron(B, eye(d));
    Cd = kron(C, eye(d));
    Dd = kron(D, eye(d));

    % Initialize storage
    x = zeros(d, num_steps+1);
    x_temp = zeros(d, num_steps+1);
    u = zeros(d, num_steps);
    % Set initial state
    x(:,1) = x0;
    x_temp(:, 1) = x0;
    % Iterate over time steps
    for k = 1:num_steps
        % Get input at time k
        u(:, k) = grad_f(Cd * x(:, k));
        % State update
        x_temp(:, k+1) = Ad * x(:, k) + Bd * u(:, k);
        y_proj = proj_fun(x_temp(:, k+1));
        x(:, k+1) = y_proj;
    end

    % Store results in struct
    results.x = x;  % State/Output trajectory
    results.x_temp = x_temp; % state before projection
    results.u = u;  % Input trajectory
end
