function dydt = model1D(t,y,paramSym)
% 1D Model equation
% paramSym: the struct of parameters

% reduced 1-D model
X = y(1); % cell population of single species

% recover S and M from X
S = paramSym.Sin - paramSym.a .* X;
M = paramSym.Min + paramSym.k .* X - paramSym.b .* X;

dXdt = (paramSym.U(S,M) - paramSym.dilution) .* X;

dydt = dXdt;

end

