function ParamSym = makeSym(Param)
% keep only base names

% Inpuy
ParamSym.Sin = Param.S1in;
ParamSym.Min = Param.M1in;

% Dilution rate (D)
ParamSym.dilution = Param.dilution;

% ratio of substrate uptake to growth (a = a1 = a2)
ParamSym.a = Param.a1;

% factor for metabolite production (k = k1 = k2)
ParamSym.k = Param.k1;

% factor for metabolite uptake (b = b1 = b2)
ParamSym.b = Param.b1;


% cell growth rate (dual Monod expression mu = mu1 = mu2)

ParamSym.U = Param.U1;

