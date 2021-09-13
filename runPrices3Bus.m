clear all
close all
clc

% This file contains the code to run the numerical experiments on the 3-bus
% network in the paper by Winnicki, Ndrio and Bose.

%---------Load the 3-bus network------------------

% For experiments 1-3, load "pglib_opf_case3_lmbdTRIAL3". For experiment 4, 
% load "pglib_opf_case3_Exp4". 

% mpc = loadcase('pglib_opf_case3_lmbdTRIAL3');
mpc = loadcase('pglib_opf_case3_Exp4');  

if isempty(find(mpc.branch(:, 6) > 0))
    disp('There are no line limits in this case')
    return;
end

%--------Define some constants to interact with the matrices---------
n           = size(mpc.bus, 1);
n_branches  = size(mpc.branch(:, 1));
n_branches  = n_branches(1);

genBuses    = mpc.gen(:, 1);

PgMax       = zeros(n, 1);
PgMin       = zeros(n, 1);

QgMax       = zeros(n, 1);
QgMin       = zeros(n, 1);

PgMax(genBuses) ...
            = mpc.gen(:, 9) / mpc.baseMVA;
PgMin(genBuses) ...
            = mpc.gen(:, 10) / mpc.baseMVA;
QgMax(genBuses) ...
            = mpc.gen(:, 4) / mpc.baseMVA;
QgMin(genBuses) ...
            = mpc.gen(:, 5) / mpc.baseMVA;
       
Pd          = mpc.bus(:, 3) / mpc.baseMVA;
Qd          = mpc.bus(:, 4) / mpc.baseMVA;

Fmax        = mpc.branch(:, 6) / mpc.baseMVA;

costGen2    = zeros(n, 1);
costGen1    = zeros(n, 1);

costGen2(genBuses) ...
            = mpc.gencost(:, 5) * mpc.baseMVA;

costGen1(genBuses) ...
            = mpc.gencost(:, 6);

WMax        = mpc.bus(:, 12) .^ 2;
WMin        = mpc.bus(:, 13) .^ 2;


%%----------Alterations-----------------

%Experiment 1  Uncomment the below line to run experiment 1
Pd =[ 0.79; 0; 1.90];
Fmax = 0.24 * ones(3, 1);
WMin = [0.95; 0.98; 0.99];
WMax = [0.98; 1.01; 1.01];

%Experiment 2 Uncomment the below lines to run experiment 2
% Pd =[ 0.79; 0; 2.00];
% Fmax = 0.20 * ones(3, 1);
% WMin = [0.95; 0.98; 0.95];
% WMax = [1.05; 1.01; 1.01];
% Qd = [0.1; 0; 0]; 

% % Experiment 3 Uncomment the below lines to run experiment 3
% Pd =[0.79; 0; 2.0];
% Fmax = 0.4 * ones(3, 1);
% WMin = [1.01; 0.98; 0.99];
% WMax = [1.05; 1.01; 1.01];

% Experiment 4 Uncomment the below lines to run experiment 4
Pd =[ 1.10; 1.10; 0.95];
Fmax = 0.9 * ones(3, 1);
WMin = [.98; .99; .95];
WMax = [1.01; 1.01; 1.02];
Qd = [1; 1; 1];

define_constants
mpc.bus(:, VMIN) = sqrt(WMin);
mpc.bus(:, VMAX) = sqrt(WMax);
mpc.branch(:, RATE_A) = Fmax;
mpc.branch(:, RATE_B) = Fmax;
mpc.branch(:, RATE_C) = Fmax;
mpc.bus(:, PD) = Pd;
mpc.gen(:, PMAX) = PgMax;
mpc.bus(:, QD) = Qd;

%-----------------Run Matpower---------------

disp('Running ACOPF with Matpower')
opt = mpoption(     'VERBOSE', 1, 'OUT_ALL', 0, ...
                    'OPF_FLOW_LIM', 1, 'OPF_ALG', 0, 'OPF_VIOLATION', 1e-6 );

results = runopf(mpc, opt);
Pg_AC =  results.gen(:, PG);
Qg_AC =  results.gen(:, QG);

cost_matpower = results.f;
disp(strcat('Matpower optimal value = ', num2str(cost_matpower)))

lambda_p_matpower = results.bus(:, LAM_P);
lambda_q_matpower = results.bus(:, LAM_Q);

disp('Prices_matpower');
disp([lambda_p_matpower, lambda_q_matpower]);

MS_ac = sum(lambda_p_matpower.*Pd - lambda_p_matpower.*Pg_AC + lambda_q_matpower.*Qd - lambda_q_matpower.*Qg_AC);
disp(strcat('MS_ac=', num2str(MS_ac)));

disp(strcat('Voltage squared=', num2str( (results.bus(:, VM).^2)')));

disp(strcat('Pgen_ac=', num2str( (results.gen(:, PG))')));
disp(strcat('Qgen_ac=', num2str( (results.gen(:, QG))')));

disp('-----------------------------------')


%------------------Run SDP-------------------


disp('Running SDP prices')
priceSDP

if strcmp(cvx_status, 'Solved') ~= 1
    disp('Problems in optimization');
    return
end

cost_sdp = cvx_optval;
disp(strcat('SDP optimal value = ', num2str(cost_sdp)))

lambda_p_sdp = lambda(1:n);
lambda_q_sdp = lambda(n+1:2*n);

disp('Prices_sdp');
disp([lambda_p_sdp, lambda_q_sdp]);

disp(strcat('Pgen_sdp=', num2str(Pg')));
disp(strcat('Qgen_sdp=', num2str((Qg'))));

MS_sdp = - sum(Pinj .* lambda_p_sdp + Qinj .* lambda_q_sdp);
disp(strcat('MS_sdp=', num2str(MS_sdp)));

disp(strcat('Eigenvalues of W = [ ', num2str(eigs(W)'), ']'))
disp(strcat('Voltage magnitude squared in SDP =', num2str((diag(W)'))));

%-----------Calculate side-payments with SDP prices and AC dispatch-----------

pi_SO_sdp = lambda_p_sdp .* Pg_AC + lambda_q_sdp .* Qg_AC - costGen2 .* (Pg_AC .* Pg_AC) - costGen1 .* Pg_AC;
        
cvx_begin quiet
    variables Pg_pi(n) Qg_pi(n);
    maximize (lambda_p_sdp' * Pg_pi + lambda_q_sdp' * Qg_pi ...
                    - costGen2' * (Pg_pi .* Pg_pi) - costGen1' * Pg_pi);
    subject to
        Pg_pi >= PgMin;
        Pg_pi <= PgMax;
        Qg_pi >= QgMin;
        Qg_pi <= QgMax;
cvx_end  


disp(strcat('Pg for max profit with SDP LMPs and AC dispatch =', num2str(Pg_pi'))); 
disp(strcat('Qg for max profit with SDP LMPs and AC dispatch =', num2str(Qg_pi'))); 

pi_max_sdp = lambda_p_sdp .* Pg_pi + lambda_q_sdp .* Qg_pi ...
                - costGen2 .* (Pg_pi .* Pg_pi) - costGen1 .* Pg_pi;

disp(strcat('Max profits with SDP LMPs =', num2str(pi_max_sdp'))); 
disp(strcat('Profits with SDP LMPs and AC dispatch =', num2str(pi_SO_sdp'))); 


side_payments_sdp = (pi_max_sdp - pi_SO_sdp);
disp(strcat('Side payments in SDP =', num2str(side_payments_sdp')));


%-----------Calculate side-payments with AC prices and AC dispatch-----------

pi_SO_ac = lambda_p_matpower .* Pg_AC + lambda_q_matpower .* Qg_AC - costGen2 .* (Pg_AC .* Pg_AC) - costGen1 .* Pg_AC;
        
cvx_begin quiet
    variables Pg_pi(n) Qg_pi(n);
    maximize (lambda_p_matpower' * Pg_pi + lambda_q_matpower' * Qg_pi ...
                    - costGen2' * (Pg_pi .* Pg_pi) - costGen1' * Pg_pi);
    subject to
        Pg_pi >= PgMin;
        Pg_pi <= PgMax;
        Qg_pi >= QgMin;
        Qg_pi <= QgMax;
cvx_end  

pi_max_ac = lambda_p_matpower .* Pg_pi + lambda_q_matpower .* Qg_pi ...
                - costGen2 .* (Pg_pi .* Pg_pi) - costGen1 .* Pg_pi;

side_payments_ac = (pi_max_ac - pi_SO_ac);
disp(strcat('Side payments in AC =', num2str(side_payments_ac')));


