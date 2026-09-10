function c = solveFixedPoint(Param)

% Solve the fixed points
syms X1 X2

S1 = Param.S1in - X1*Param.a1;
S2 = Param.S2in - X2*Param.a2;
M1 = Param.M1in - Param.b2 * X2 + Param.k1 * X1;
M2 = Param.M2in - Param.b1 * X1 + Param.k2 * X2;

Fx1 = (Param.U1(S1,M2) - Param.dilution).* X1;
Fx2 = (Param.U2(S2,M1) - Param.dilution).* X2;
jacob = jacobian([Fx1; Fx2], [X1, X2]);

jacobFun = matlabFunction(jacob, 'Vars',[X1, X2]);

sol = solve([Fx1 == 0; Fx2 == 0], [X1 X2]);

solution = double([sol.X1, sol.X2]);
numSol = size(solution,1);

type = NaN * ones(numSol,1); % NaN: invalid/infeasible solution; 1: physically meaningful solution

for i = 1: numSol
    x0 = solution(i,1);
    y0 = solution(i,2);
    solS1 = Param.S1in - x0*Param.a1;
    solS2 = Param.S2in - y0*Param.a2;
    solM1 = Param.M1in - Param.b2 * y0 + Param.k1 * x0;
    solM2 = Param.M2in - Param.b1 * x0 + Param.k2 * y0;
    
    solT = [x0, y0, solS1, solS2, solM1, solM2];
    if any(solT<0 | imag(solT) ~= 0)
        continue;
    end

    type(i) = 1;
end

effSol = solution(type == 1, :);

numEff = size(effSol, 1);

% Determine the stability of each fixed point and calculate separatrix
c = cell(numEff, 3);
    % three columns:
    % Col 1: solution
    % Col 2: type (0: unstable; 1: stable; 2:non-hyperbolic)
    % Col 3: separatrix for saddle points

tolE = 0;
tolF = 0;
tolB = 0;

for i = 1:numEff
    x0 = effSol(i,1);
    y0 = effSol(i,2);
    c{i,1} = effSol(i,:);

    J = jacobFun(x0, y0);
    [V,D] = eig(J);
    eigVal = real(diag(D));

    if any(eigVal>tolE)
        c{i,2} = 0; % unstable

        if any(eigVal<-tolE) 
            % compute separatrix for the saddle point

            [~, idx] = min(eigVal);
            v_stable = V(:,idx) / norm(V(:,idx)); 
            delta = 1e-4;
            z_init1 = [x0; y0] + delta * v_stable;
            z_init2 = [x0; y0] - delta * v_stable;
            
            % backward
            odefun = @(t, z) -[(Param.U1(Param.S1in - z(1)*Param.a1, Param.M2in - Param.b1 * z(1) + Param.k2 * z(2)) - Param.dilution) .* z(1); ...
                               (Param.U2(Param.S2in - z(2)*Param.a2, Param.M1in - Param.b2 * z(2) + Param.k1 * z(1)) - Param.dilution) .* z(2)];  
            
            
            opt = odeset('Events', @(t,z) stopAtBoundary(t,z,Param,tolB));
    
            if isFeasible(z_init1, Param, tolF)
                [~, Z1] = ode45(odefun, 0:0.01:50, z_init1, opt);
            else
                Z1 = [];
            end
            if isFeasible(z_init2, Param, tolF)
                [~, Z2] = ode45(odefun, 0:0.01:50, z_init2, opt);
            else
                Z2 = [];
            end
            Zt = [Z1(end:-1:1,:);Z2];
    
            c{i,3} = Zt;
        end
        
    elseif all(eigVal < -tolE)
        c{i,2} = 1; % stable
        
    else
        c{i,2} = 2; % non-hyperbolic / near-critical
    end
end

end


function tf = isFeasible(z, Param, tol)

X1 = z(1);
X2 = z(2);

S1 = Param.S1in - Param.a1*X1;
S2 = Param.S2in - Param.a2*X2;
M1 = Param.M1in - Param.b2*X2 + Param.k1*X1;
M2 = Param.M2in - Param.b1*X1 + Param.k2*X2;

tf = all([X1, X2, S1, S2, M1, M2] >= -tol);

end

function [value, isterminal, direction] = stopAtBoundary(~, z, Param, tol)

X1 = z(1);
X2 = z(2);

S1 = Param.S1in - Param.a1*X1;
S2 = Param.S2in - Param.a2*X2;
M1 = Param.M1in - Param.b2*X2 + Param.k1*X1;
M2 = Param.M2in - Param.b1*X1 + Param.k2*X2;

value = [X1; X2; S1; S2; M1; M2] + tol;
isterminal = ones(6,1);
direction = -ones(6,1);

end