function results = projected_triple_momentum_2(A, B, C, D, x0, num_steps, grad_f, proj_fun, P)
    n = length(x0);
    d = n/2;
    % Kronecker
    Ad = kron(A, eye(d));
    Bd = kron(B, eye(d));
    Cd = kron(C, eye(d));
    Dd = kron(D, eye(d));

    % Get system dimensions
    n = size(Ad, 1); % State dimension
    m = size(Bd, 2); % Input dimension
    p = size(Cd, 1); % Output dimension

     % Initialize storage
    x = zeros(n, num_steps+1);
    x_temp = zeros(n, num_steps+1);
    y = zeros(p, num_steps);
    u = zeros(m, num_steps);

    % Set initial state
    x(:,1) = x0;
    x_temp(:,1) = x0;

    for k = 1:num_steps
        u(:, k) = grad_f(Cd * x(:, k));
        x_temp(:, k+1) = Ad * x(:, k) + Bd * u(:, k);
        y_proj = proj_fun(x_temp(1:d, k+1));
        x(1:d, k+1) = y_proj;
        Shift = kron(P(2:end,2:end)\P(1,2:end)', eye(d))*(y_proj - x_temp(1:d, k+1));
        x(d+1:end, k+1) = x_temp(d+1:end, k+1) - Shift(1:d);
        y(:, k) = y_proj;
    end

    results.x = x;
    results.y = y;
    results.u = u;
end
