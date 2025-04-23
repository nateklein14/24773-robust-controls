function [inf_norm] = calc_inf_norm(sys, tol)

gamma_h = 10000;
gamma_l = 0;
gamma_new = (gamma_h + gamma_l)/2;

G = minreal(ss(sys));
dim = size(G.D, 1);
A = G.A;
B = G.B;
C = G.C;
D = G.D;

while(abs(gamma_h - gamma_l) > tol)
    gamma_try = gamma_new;
    if max(svds(D)) >= gamma_try % initial check
        gamma_l = gamma_try;
    else
        R = eye(dim)*gamma_try^2 - D'*D;
        Abar = A + B/R*D'*C;
        H = [Abar B/R*B' ; -C'*(eye(dim) + D/R*D')*C -Abar']; % hamiltonian
        e = eig(H);
        if any(abs(real(e)) < 1e-14) % H has an imaginary eigenvalue
            gamma_l = gamma_try; % norm is not less than gamma
        else
            gamma_h = gamma_try; % norm is less than gamma
        end
    end
    gamma_new = (gamma_h + gamma_l)/2;
end

inf_norm = gamma_new;
end