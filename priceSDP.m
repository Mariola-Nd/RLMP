% This computes the SDP prices, given the matpower file with the case.

[Ybus, ~, ~] ...
            = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);

        
e_mat = eye(n);

Phi         = {};
Psi         = {};
JJ          = {};

j           = sqrt(-1);

for k=1:n
    y_k         = e_mat(:,k) * (e_mat(:, k)') * Ybus;
    Phi{k}      = (ctranspose(y_k) + y_k)/2;
    Psi{k}      = (ctranspose(y_k) - y_k)/(2*j);
end


clear y_k

Yf = {};
Yt = {};

for bb = 1:n_branches

    Yf{bb}  = zeros(n);
    Yt{bb}  = zeros(n);

    ii = mpc.branch(bb, 1);
    jj = mpc.branch(bb, 2);
    
    Yf{bb}(jj, ii) ...
            = conj(Ybus(ii, jj));
    Yf{bb}(ii, ii) ...
            = conj(-Ybus(ii, jj));
    Yf{bb}  = (Yf{bb} + ctranspose(Yf{bb}))/2;

    jj = mpc.branch(bb, 1);
    ii = mpc.branch(bb, 2);
    
    Yt{bb}(jj, ii) ...
            = conj(Ybus(ii, jj));
    Yt{bb}(ii, ii) ...
            = conj(-Ybus(ii, jj));
    Yt{bb}  = (Yt{bb} + ctranspose(Yt{bb}))/2;
    
end

cvx_solver sedumi
cvx_precision high

cvx_begin quiet
    variables Pg(n) Qg(n) Pinj(n) Qinj(n); 
    variables Pf(n_branches)  Pt(n_branches) aux(n);
    variable W(n, n) hermitian;
    dual variable lambda;
    dual variables mu_kl{2 * n_branches};
    dual variable U;
    minimize sum(aux);
    subject to
                
        for kk = 1:n
            Pinj(kk) == real( trace( Phi{kk} * W ));
            Qinj(kk) == real( trace( Psi{kk} * W ));
        end
    
        lambda : [Pg - Pd - Pinj; Qg - Qd - Qinj] == zeros(2*n, 1);
        
        Pg <= PgMax;
        Pg >= PgMin;
        Qg <= QgMax;
        Qg >= QgMin;
       
        real(diag(W)) >= WMin;
        real(diag(W)) <= WMax;
        
        % Cost construction:
        for cc = 1:n
            costGen2(cc) * square(Pg(cc)) + costGen1(cc) * Pg(cc) <= aux(cc);
        end
            
        % Line limits
        for bb = 1:n_branches
            mu_kl{bb} : real(trace(Yf{bb} * W)) <= Fmax(bb);
            mu_kl{bb + n_branches} : real(trace(Yt{bb} * W)) <= Fmax(bb);    
        end
        
        U : W == hermitian_semidefinite( n );

cvx_end

