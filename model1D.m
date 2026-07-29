function dydt = model1D(t,y,paramSym)
% 1D Model equation
% paramSym: the struct of parameters
% inputSym: the struct of input substrate/metabolite concentrations

% reduced 1-D model
X = y(1); % cell population of single species

% recover S1,S2,M1,M2 from X1,X2
S = paramSym.Sin - paramSym.a .* X;
M = paramSym.Min + paramSym.k .* X - paramSym.b .* X;

dXdt = (paramSym.U(S,M) - paramSym.dilution) .* X;

dydt = dXdt;

end

