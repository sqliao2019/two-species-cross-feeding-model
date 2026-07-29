mumax_perturbed = 4.2;
% Normal:                   5.0;
% Weak perturbation:        4.2;
% Moderate perturbation:    3.0;
% Strong perturbation:      1.0;

Param_normal = makeSym(getPara(5));
Param_perturb = makeSym(getPara(mumax_perturbed));

%%
opt = odeset('NonNegative',1);

T0 = 0; 
T1 = 3;
T2 = 4;
T3 = 20;
T4 = 20;

X_ini = 2;

y0 = X_ini;

[tb,yb] = ode45(@(t,y)model1D(t,y,Param_normal), [T0 T1], y0,       opt);

[td,yd] = ode45(@(t,y)model1D(t,y,Param_perturb), [T1 T2], yb(end,:),opt);
[ta,ya] = ode45(@(t,y)model1D(t,y,Param_normal), [T2 T3], yd(end,:),opt);

[tp,yp] = ode45(@(t,y)model1D(t,y,Param_perturb), [T1 T4], yb(end,:),opt);

figure;hold on;
plot(tb,yb(:,1));
plot(td,yd(:,1));
plot(ta,ya(:,1));

figure;hold on;
plot(tb,yb(:,1));
plot(tp,yp(:,1));

%% 1-D Landscape (normal)

figure;hold on;

Xmax = Param_normal.Sin/Param_normal.a;

[X1D,U1D] = ode45(@(X,U)Gradient(X,U,Param_normal), linspace(0,Xmax,1000), 0);
U1D = (U1D-U1D(1))./(max(U1D)-min(U1D));
plot(X1D, U1D);

markerX_1 = yb(1,1);    
markerY_1 = interp1(X1D, U1D, markerX_1);

markerX_2 = yb(end,1);
markerY_2 = interp1(X1D, U1D, markerX_2);

markerX_3 = yd(end,1);
markerY_3 = interp1(X1D, U1D, markerX_3);

markerX_4 = ya(end,1);
markerY_4 = interp1(X1D, U1D, markerX_4);

plot(markerX_1, markerY_1, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','k');
plot(markerX_2, markerY_2, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','r');
plot(markerX_3, markerY_3, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','b');
plot(markerX_4, markerY_4, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','g');
xlim([0 Xmax])
ylim([min(U1D)-0.2 max(U1D)+0.2])
box on;
yticks([]);
xticks([0 Xmax])
%xticklabels(["0", "+"])


%% 1-D Landscape (perturbed)

figure;hold on;

Xmax = Param_perturb.Sin/Param_perturb.a;

[X1D,U1D] = ode45(@(X,U)Gradient(X,U,Param_perturb), linspace(0,Xmax,1000), 0);
U1D = (U1D-U1D(1))./(max(U1D)-min(U1D));
plot(X1D, U1D);

markerX_2 = yd(1,1);    
markerY_2 = interp1(X1D, U1D, markerX_2);

markerX_3 = yd(end,1);
markerY_3 = interp1(X1D, U1D, markerX_3);

markerX_5 = yp(end,1);
markerY_5 = interp1(X1D, U1D, markerX_5);

plot(markerX_2, markerY_2, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','r');
plot(markerX_3, markerY_3, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','b');
plot(markerX_5, markerY_5, 'o', 'MarkerEdgeColor','none','MarkerFaceColor','m');

xlim([0 Xmax])
ylim([min(U1D)-0.2 max(U1D)+0.2])
box on;
yticks([]);
xticks([0 Xmax])
% xticklabels(["0", "+"])
