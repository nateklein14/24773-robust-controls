function [K, CL, GAM] = myhinfsyn(P, nmeas, ncont)

sys = ss(P);
A = sys.A;
B1 = sys.B(:,1:size(sys.B,2)-ncont);
B2 = sys.B(:,size(sys.B,2)-ncont+1:end);
C1 = sys.C(1:size(sys.C,1)-nmeas,:);
C2 = sys.C(size(sys.C,1)-nmeas+1:end,:);
D12 = sys.D(1:size(sys.D,1)-ncont, size(sys.D,2)-nmeas+1:end);
D21 = sys.D(size(sys.D,1)-ncont+1:end, 1:size(sys.D,2)-nmeas);

gamma_h = 10000;
gamma_l = 0;
gamma_try = (gamma_h+gamma_l)/2;

while abs(gamma_h-gamma_l) > 1e-6
    flag = 0;
    gamma_new = gamma_try;
    Aare = A;
    Bare = B2*inv(D12'*D12)*B2' - gamma_try^-2*B1*B1';
    Care = C1'*C1;
    try
        X = are(Aare, Bare, Care);
        try
            Y = are(Aare', C2'*inv(D21*D21')*C2 - gamma_try^-2*C1'*C1, B1*B1');
            if max(abs(eig(X*Y))) < gamma_try^2
                gamma_h = gamma_try;
            else
                flag = 1;
            end
        catch
            flag = 1;
        end
    catch
        flag = 1;
    end
    if flag
        gamma_l = gamma_try;
    end
    gamma_try = (gamma_h + gamma_l)/2;
end

GAM = gamma_try;
Linf = Y*inv(eye(size(X*Y))-GAM^-2*X*Y)*(C2+GAM^-2*D21*B1'*X)'*inv(D21*D21');
K = ss(A+GAM^-2*B1*B1'*X-B2*inv(D12'*D12)*B2'*X-Linf*(C2+GAM^-2*D21*B1'*X), ...
       Linf, ...
       -inv(D12'*D12)*B2'*X, ...
       0);
CL = lft(P, K);

end