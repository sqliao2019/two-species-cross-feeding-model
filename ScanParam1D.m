function Sol = ScanParam1D(param,scanName,varRange, arrayLen)

syms x v

scanName = char(string(scanName));

param.(scanName) = v;

Xsol = solve(Gradient(x,0,param)==0, x,"ReturnConditions",true);

XsolFun = matlabFunction(Xsol.x,'var',v);

% Jacobian of dX/dt; in one dimension, this is also the eigenvalue
jacobianFun = matlabFunction( ...
    -diff(Gradient(x,0,param), x), ...
    'Vars', [x v]);

varScan = linspace(varRange(1), varRange(2), arrayLen);

numBranches = numel(Xsol.x);
Sol = NaN(numBranches, arrayLen);

Lambda = NaN(size(Sol));

for i = 1:arrayLen
    sol_test = XsolFun(varScan(i));
    param.(scanName) = varScan(i);
    for j = 1:length(sol_test)
        X_test = sol_test(j);

        S1 = param.Sin - param.a*X_test;
        M2 = param.Min - param.b*X_test + param.k*X_test;
        if (~isreal(X_test) || S1<0 || M2<0 || X_test<0)
            sol_test(j) = NaN;
            continue;
        end
        Lambda(j,i) = jacobianFun(X_test, varScan(i));
    end

    Sol(:,i) = sol_test;
end

Sol(abs(imag(Sol))>0)=NaN;

% Classify each equilibrium point according to its eigenvalue
stableMask   = isfinite(Sol) & isfinite(Lambda) & Lambda < 0;
unstableMask = isfinite(Sol) & isfinite(Lambda) & Lambda > 0;
neutralMask  = isfinite(Sol) & isfinite(Lambda) & Lambda == 0;

StableSol = NaN(size(Sol));
UnstableSol = NaN(size(Sol));
NeutralSol = NaN(size(Sol));

StableSol(stableMask) = Sol(stableMask);
UnstableSol(unstableMask) = Sol(unstableMask);
NeutralSol(neutralMask) = Sol(neutralMask);

stableColor   = [0 0.4470 0.7410];
unstableColor = [0.8500 0.3250 0.0980];
mixedColor    = [0.5 0.5 0.5];

figure;
hold on;

for j = 1:size(Sol,1)
    plot(varScan, StableSol(j,:), ...
        'Color', stableColor, ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    plot(varScan, UnstableSol(j,:), ...
        'Color', unstableColor, ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    plot(varScan, NeutralSol(j,:), ...
        'o', ...
        'Color', mixedColor, ...
        'MarkerFaceColor', mixedColor, ...
        'LineStyle', 'none');
end

xlim(varRange)
ylim([-1 6])


end

