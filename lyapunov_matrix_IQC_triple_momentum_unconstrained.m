function [P,K,delta,epsi] = lyapunov_matrix_IQC_triple_momentum_unconstrained(A, B, C, D, m, L)
% yalmip('version') 
% ans =
% 
%     '20230622'
% mosekopt('version')     
% MOSEK Version 10.2.11 (Build date: 2024-12-3 13:30:04)
% Copyright (c) MOSEK ApS, Denmark WWW: mosek.com
% Platform: Windows/64-X86

%%
yalmip('clear')
% parameters
rho = 1 - 1/(sqrt(L/m)); % fastest convergence rate
%% algorithm analysis using IQC
% sector IQC
Dphiy1 = [L; -m];
Dphiu1 = [-1; 1];

% off-by-one IQC
Aphi2 = [0];
Bphiy2 = [-L];
Bphiu2 = [1];
Cphi2 = [1; 0];
Dphiy2 = [L; -m];
Dphiu2 = [-1; 1];

% weighted off-by-one IQC
rho_bar = rho;
Aphi3 = [0];
Bphiy3 = [-L];
Bphiu3 = [1];
Cphi3 = [rho_bar^2; 0];
Dphiy3 = [L; -m];
Dphiu3 = [-1; 1];

M = [0 1; 1 0];
%% Augmented system matrices
Ahat = [A zeros(2,1) zeros(2,1);...
        Bphiy2*C Aphi2 0;...
        Bphiy3*C 0 Aphi3];
%
Bhat = [B;...
        Bphiy2*D+Bphiu2;...
        Bphiy3*D+Bphiu3];
%
Chat = [Dphiy1*C zeros(2,1) zeros(2,1);...
        Dphiy2*C Cphi2 zeros(2,1);...
        Dphiy3*C zeros(2,1) Cphi3];
%
Dhat = [Dphiy1*D+Dphiu1;...
        Dphiy2*D+Dphiu2;...
        Dphiy3*D+Dphiu3];
%%
% sdp variables
a1 = sdpvar;
a2 = sdpvar;
a3 = sdpvar;
P = sdpvar(4,4);
epsi = sdpvar;
delta = sdpvar;
K = [Ahat'*P*Ahat-rho^2*P, Ahat'*P*Bhat; Bhat'*P*Ahat, Bhat'*P*Bhat] + [Chat Dhat]'*kron(diag([a1, a2, a3]),M)*[Chat Dhat];
Constraints = [P>=delta*eye(4); K<=-epsi*eye(5); a1>=0; a2>=0; a3>=0; epsi>=0; delta>=0];
ops = sdpsettings('solver','mosek','debug',1,'verbose',0);
ops.mosek.MSK_DPAR_INTPNT_TOL_REL_GAP    = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_PFEAS      = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_DFEAS      = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_MU_RED     = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_PATH       = 1e-10;
ops.mosek.MSK_DPAR_INTPNT_TOL_STEP_SIZE  = 1e-10;
sol = optimize(Constraints,[],ops);

if sol.problem ~= 0
    warning(['Lyapunov matrix searching status: ', sol.info]);
    % yalmiperror(sol.problem)
end
% 
P = value(P);
K = value(K);
epsi = value(epsi);
delta = value(delta);
end