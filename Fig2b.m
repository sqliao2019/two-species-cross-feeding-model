Param = getPara(5);

numVF = 16;
x1v = linspace(0, Param.S1in/Param.a1, numVF);
x2v = linspace(0, Param.S2in/Param.a2, numVF);

[X1,X2] = meshgrid(x1v, x2v);

S1 = Param.S1in - Param.a1 .* X1;
S2 = Param.S2in - Param.a2 .* X2;
M1 = Param.M1in + Param.k1 .* X1 - Param.b2 .* X2;
M2 = Param.M2in + Param.k2 .* X2 - Param.b1 .* X1;
mask = S1<-1e-5 | S2<-1e-5 | M1<-1e-5 | M2<-1e-5;

Fx1 = (Param.U1(S1,M2) - Param.dilution).* X1;
Fx2 = (Param.U2(S2,M1) - Param.dilution).* X2;

Fx1(mask) = NaN;
Fx2(mask) = NaN;

mag = sqrt(Fx1.^2 + Fx2.^2);

Fx1 = Fx1./mag;
Fx2 = Fx2./mag;

magNorm = (mag - min(mag(:))) ./ (max(mag(:)) - min(mag(:)));
cmap = parula(256);

figure;
axis equal; hold on;

for i = 1:numel(X1)
    % Pick color based on magnitude
    idx = max(1, round(magNorm(i)*255)+1);
    thisColor = cmap(idx,:);
    
    % Plot each arrow with that color
    h = quiver(X1(i), X2(i), Fx1(i), Fx2(i), 0.15, ...
           'Color', thisColor, 'MaxHeadSize', 3,'LineWidth', 1.2);
    h.ShowArrowHead = 'on'; 

end

edgewidth = 0.00;
axis([0-edgewidth, Param.S1in/Param.a1+edgewidth, 0-edgewidth, Param.S2in/Param.a2+edgewidth]);

colormap(cmap);
cb = colorbar;
cb.Label.String = 'Vector magnitude';

plot(x1v, Param.b1 ./Param.k2 .*x1v, 'k');
plot(Param.b2 ./Param.k1 .*x1v, x2v, 'k');

fixedPoint = solveFixedPoint(Param);
plot(fixedPoint{1,1}(1), fixedPoint{1,1}(2), 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[0.1 0.1 0.1]);
plot(fixedPoint{2,3}(:,1), fixedPoint{2,3}(:,2),'k');
plot(fixedPoint{2,1}(1), fixedPoint{2,1}(2), 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[1 1 1]);
plot(fixedPoint{3,1}(1), fixedPoint{3,1}(2), 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[0.1 0.1 0.1]);