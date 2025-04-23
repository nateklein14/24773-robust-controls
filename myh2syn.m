function [K, CL, GAM] = myh2syn(P, nmeas, ncont)

sys = ss(P);
A = sys.A;
B1 = sys.B(:,1:size(sys.B,2)-ncont);
B2 = sys.B(:,size(sys.B,2)-ncont+1:end);
C1 = sys.C(1:size(sys.C,1)-nmeas,:);
C2 = sys.C(size(sys.C,1)-nmeas+1:end,:);
D12 = sys.D(1:size(sys.D,1)-ncont, size(sys.D,2)-nmeas+1:end);
D21 = sys.D(size(sys.D,1)-ncont+1:end, 1:size(sys.D,2)-nmeas);

X = are(A-B2*D12'*C1, B2*B2', C1'*(eye(size(D12,1))-D12*D12')*C1);
Y = are((A-B1*D21'*C2)', C2'*C2, B1*(eye(size(D21,2))-D21'*D21)*B1');
F = -B2'*X+D12'*C1;
L = -Y*C2'+B1*D21';

K = ss(A+B2*F+L*C2, -L, F, 0);
GAM = sqrt(trace(B1'*X*B1) + trace(F*Y*F'));
CL = lft(P, K);

end