%% Convergence performance comparison among projected gradient, projected triple momentum Algorithm 1 and 2
% ver
% -----------------------------------------------------------------------------------------------------
% MATLAB Version: 24.1.0.2689473 (R2024a) Update 6
% Operating System: Microsoft Windows 11 Pro Version 10.0 (Build 26100)
% Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% -----------------------------------------------------------------------------------------------------

% yalmip('version') 
% ans =
% 
%     '20230622'

% mosekopt('version')     
% MOSEK Version 10.2.11 (Build date: 2024-12-3 13:30:04)
% Copyright (c) MOSEK ApS, Denmark WWW: mosek.com
% Platform: Windows/64-X86

clear all;
%% system matrices
syms alpha beta gamma positive;
% Triple Momentum algorithm original state-space representation
A_O = [1+beta -beta; 1 0];
B_O = [-alpha; 0];
C_O = [1+gamma -gamma];
% linear transformation that does not change the algorithm
T = [1+gamma -gamma; 0 1];
V = T^(-1);
A_T = simplify(T*A_O*V);
B_T = T*B_O;
C_T = C_O*V;
D_T = zeros(size(C_T,1), size(B_T,2));
%% Objective function and parameters
% Define objective function
d = 2;
F = [100 -1; -1 1];
p = [1; 10];
f_handle = @(x) 0.5 * x' * F * x + p' * x;
% d = 3;
% F =  [100 -1 1; -1 1 2; 1 2 10];
% p = [1; 10; 5];
% f_handle = @(x) 0.5 * x' * F * x + p' * x;

% Compute gradient, strong convexity constant and Lipchitz constant (for
% quadratic function only) 
% strong convexity constant m = 0.9899;  and
% Lipschitz constant L = 100.0101;   in our example.
[grad_f, m, L] = analyze_quadratic_function(f_handle, d);
%% Algorithm system matrices
rho = 1 - 1/(sqrt(L/m)); % fastest convergence rate
% parameters of the TM method
alpha_opt = (1 + rho)/L;
beta_opt = rho^2/(2-rho);
gamma_opt = rho^2/((1+rho)*(2-rho));
% numerical state-space matrices of the TM method given m and L
A = double(subs(A_T, [alpha, beta, gamma], [alpha_opt, beta_opt, gamma_opt]));
B = double(subs(B_T, [alpha, beta, gamma], [alpha_opt, beta_opt, gamma_opt]));
C = double(C_T);
D = zeros(size(C_T,1), size(B_T,2));
%% Constraints
% Ellipsoid constraint
Q = [1 0; 0 2];
% all constraints
constraint_fun = @(y) { y' * Q * y <= 1;
                        % y(1) >= 0;
                        % y(2) <= 0.5
                        };

% Q = [1 0 0; 0 2 0; 0 0 3];
% constraint_fun = @(y) { 
%     y' * Q * y <= 1;
%                         % y(1) >= 5;
%                         % y(2) >= -1;
%                         % y(3) >= 0
%                         };

% Define projection step
proj_fun = build_projection(constraint_fun);  % or with other constraints
%% RUN algorithms
x0 = 10*rand(2*d,1);
num_steps = 100;
% Run algorithm 1
results_proj_1 = projected_triple_momentum_1(A, B, C, D, x0, num_steps, grad_f, proj_fun);
% Compute Lyapunov function matrix P
P = lyapunov_matrix_IQC_triple_momentum_unconstrained(A, B, C, D, m, L);
% Run algorithm 2
results_proj_2 = projected_triple_momentum_2(A, B, C, D, x0, num_steps, grad_f, proj_fun, P);
%
%% Compute the optimal solution using yalmip
x_proj = sdpvar(d,1);
Constraints = [];
c_list = constraint_fun(x_proj);
Constraints = [Constraints; c_list{:}];
Objective = f_handle(x_proj);
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
ops.mosek.MSK_DPAR_INTPNT_TOL_REL_GAP    = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_PFEAS      = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_DFEAS      = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_MU_RED     = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_PATH       = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_STEP_SIZE  = 1e-10;
% ops = sdpsettings('solver', 'sedumi','sedumi.eps', 1e-12);

sol = optimize(Constraints,Objective, ops);
x_optimal_proj = value(x_proj);
f_optimal_proj = f_handle(x_optimal_proj);
%% projected gradient algorithm (to show the solver precision limit)
A_g = 1;
B_g = - 2/(L + m);
C_g = 1;
D_g = 0;
x0_g = x0(1:d);
results_proj_g = proj_gradient(A_g, B_g, C_g, D_g, x0_g, num_steps, grad_f, proj_fun);
%% Plot variable norm
figure;
semilogy(vecnorm(results_proj_g.x-x_optimal_proj, 2, 1),'LineWidth',2); hold on;
semilogy(vecnorm(results_proj_1.y-x_optimal_proj, 2, 1),'LineWidth',2); hold on;
semilogy(vecnorm(results_proj_2.y-x_optimal_proj, 2, 1),'LineWidth',2);
grid on;
% xlim([0, 40]);
% ylim([1e-13, 1e2]);
% yticks([1e-13 1e-10 1e-5 1e0]);
legend('Projected gradient descent','Algorithm 1', 'Algorithm 2','fontsize', 11);
xlabel ('Iteration Step','fontsize', 15);
ylabel('$\left\| y_{k} - y_{\Omega}^{\textup{opt}} \right\|_2$', 'fontsize', 15, 'Interpreter','latex');

%% plot objective function error
f_record_g = zeros(1, num_steps);
f_record_1 = zeros(1, num_steps);
f_record_2 = zeros(1, num_steps);
for k = 1:num_steps
    f_record_g(k) = f_handle(results_proj_g.x(:,k));
    f_record_1(k) = f_handle(results_proj_1.y(:,k));
    f_record_2(k) = f_handle(results_proj_2.y(:,k));
end
figure;
semilogy(vecnorm(f_record_g-f_optimal_proj, 2, 1),'LineWidth',2); hold on;
semilogy(vecnorm(f_record_1-f_optimal_proj, 2, 1),'LineWidth',2); hold on;
semilogy(vecnorm(f_record_2-f_optimal_proj, 2, 1),'LineWidth',2);
grid on;
% xlim([0, 40]);
% ylim([1e-13, 1e3]);
% yticks([1e-13 1e-10 1e-5 1e0]);
legend('Projected gradient descent','Algorithm 1', 'Algorithm 2','fontsize', 11);
xlabel ('Iteration Step','fontsize', 15);
ylabel('$\left| f(y_{k}) - f(y_{\Omega}^{\textup{opt}}) \right|$', 'fontsize', 15, 'Interpreter','latex');
