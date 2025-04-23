function [two_norm] = calc_2_norm(sys)

G = minreal(ss(sys));
if isstable(G) && all(G.D == 0, 'all')
    Q = gram(G, 'o');
    two_norm = sqrt(trace(G.B'*Q*G.B));
else
    two_norm = Inf;
end

end