function [grad_f, m, L] = analyze_quadratic_function(f_handle, d)
    if d == 1
        x = sym('x', 'real');
    else
        x = sym('x', [d 1], 'real');
    end
    f_sym = f_handle(x);
    grad_sym = gradient(f_sym, x);
    H = hessian(f_sym, x);
    eigenvalues = eig(H);
    m = real(double(min(eigenvalues)));  % strong convexity
    L = real(double(max(eigenvalues)));  % Lipschitz gradient
    grad_f = matlabFunction(grad_sym, 'Vars', {x});


end
