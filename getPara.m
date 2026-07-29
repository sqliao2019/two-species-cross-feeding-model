function Param = getPara(x)
%   PARAM = GETPARA(X) returns a struct of parameters and derived function 
%   handles for the model. X provides the numeric value(s) for internal
%   parameters that are used inside anonymous expressions.

% Input
Param.S1in = 5;
Param.S2in = 5; 
Param.M1in = 0;
Param.M2in = 0;

% Dilution rate (D)
Param.dilution = 1;

% ratio of substrate uptake to growth (a1, a2)
Param.a1 = 1.0;
Param.a2 = 1.0;

% factor for metabolite production (k1, k2)
Param.k1 = 0.25; 
Param.k2 = Param.k1;

% factor for metabolite uptake (b1, b2)
Param.b1 = 0.05;
Param.b2 = Param.b1;


% cell growth rate (dual Monod expression mu1, mu2)
U1max = x; 
KS1 = 1; 
KM2 = 1;
Param.U1 = @(S1, M2) U1max .* S1./(S1 + KS1) .* M2./(M2 + KM2);

U2max = x; 
KS2  = KS1; 
KM1  = KM2;
Param.U2 = @(S2, M1) U2max .* S2./(S2 + KS2) .* M1./(M1 + KM1);
