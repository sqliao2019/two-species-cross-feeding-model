function dUdX = Gradient(X, U, paramSym)
% Used for 1D gradient integral
% require assumption of symmetry

% recover S,M from X
S = paramSym.Sin - paramSym.a .* X;
M = paramSym.Min + paramSym.k .* X - paramSym.b .* X;

dUdX = - (paramSym.U(S,M) - paramSym.dilution) .* X;

end