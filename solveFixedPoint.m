function c = solveFixedPoint(Param)


% Param = getPara();

syms X1 X2

S1 = Param.S1in - X1*Param.a1;
S2 = Param.S2in - X2*Param.a2;
M1 = Param.M1in - Param.b2 .* X2 + Param.k1 .* X1;
M2 = Param.M2in - Param.b1 .* X1 + Param.k2 .* X2;

Fx1 = (Param.U1(S1,M2) - Param.dilution).* X1;
Fx2 = (Param.U2(S2,M1) - Param.dilution).* X2;
jacob = jacobian([Fx1; Fx2], [X1, X2]);

jacobFun = matlabFunction(jacob, 'Vars',[X1, X2]);

sol = solve([Fx1 == 0; Fx2 == 0], [X1 X2]);

solution = double([sol.X1, sol.X2]);
numSol = size(solution,1);

type = NaN * ones(numSol,1); % NaN: not real solution; 1: meaningful solution;

for i = 1: numSol
    x0 = solution(i,1);
    y0 = solution(i,2);
    solS1 = Param.S1in - x0*Param.a1;
    solS2 = Param.S2in - y0*Param.a2;
    solM1 = Param.M1in - Param.b2 .* y0 + Param.k1 .* x0;
    solM2 = Param.M2in - Param.b1 .* x0 + Param.k2 .* y0;
    
    solT = [x0, y0, solS1, solS2, solM1, solM2];
    if any(solT<0 | imag(solT) ~= 0)
        continue;
    end

    type(i) = 1;
end

effSol = solution(type == 1, :);

numEff = size(effSol, 1);

c = cell(numEff, 3);
% three columns:
% Col 1: solution
% Col 2: type (0: (0,0); 1: unstable; 2: stable)
% Col 3: separatrix for unstable one (assumed to be the saddle point)

for i = 1:numEff
    x0 = effSol(i,1);
    y0 = effSol(i,2);
    c{i,1} = effSol(i,:);
    if x0 == 0 && y0 == 0
        c{i,2} = 0;
    else
        J = jacobFun(x0, y0);
        [V,D] = eig(J);
        if any(real(diag(D))>0)
            c{i,2} = 1;
            [~, idx] = min(real(diag(D)));
            v_stable = V(:,idx) / norm(V(:,idx)); 
            delta = 1e-4;
            z_init1 = [x0; y0] + delta * v_stable;
            z_init2 = [x0; y0] - delta * v_stable;
            
            % backward
            odefun = @(t, z) -[(Param.U1(Param.S1in - z(1)*Param.a1, Param.M2in - Param.b1 .* z(1) + Param.k2 .* z(2)) - Param.dilution) .* z(1); ...
                               (Param.U2(Param.S2in - z(2)*Param.a2, Param.M1in - Param.b2 .* z(2) + Param.k1 .* z(1)) - Param.dilution) .* z(2)];  
           
            [~, Z1] = ode45(odefun, 0:0.01:30, z_init1);
            [~, Z2] = ode45(odefun, 0:0.01:30, z_init2);
            Zt = [Z1(end:-1:1,:);Z2];
            mask = (Param.M1in - Param.b2 .* Zt(:,2)) + Param.k1 .* Zt(:,1)>=-1e-3 &...
                   (Param.M2in - Param.b1 .* Zt(:,1) + Param.k2 .* Zt(:,2))>=-1e-3 ;
            Zt = Zt(mask,:);
            c{i,3} = Zt;
            
        else
            c{i,2} = 2;
        end
    end
end