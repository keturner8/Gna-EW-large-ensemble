c1 = load('terr_r45l18t08d05_ctrl.mat');
w1f = load('terr_r45l18t08d05_ewf.mat');
w1s = load('terr_r45l18t08d05_ews.mat');
%% behaviour of pCO2, T, and pH with terrestrial system
c = get(gca, 'ColorOrder');
figure;
subplot(3,2,1)
plot(c1.time, c1.chi*1e6, 'k'); hold on
plot(w1f.time, w1f.chi*1e6, 'Color', c(1,:))
plot(w1s.time, w1s.chi*1e6, '--', 'Color', c(1,:))
plot([2200 2200], [350 600], 'k--');
xlim([2000 2300])
title('(a) Atmospheric CO_2 concentrations')
ylabel('CO_2 (ppm)')
legend('Control, RCP 4.5, high \lambda', ...
    '2 PgC/yr EW \times 100 yr', ...
    '1 PgC/yr EW \times 200 yr')


subplot(3,2,2)
plot(c1.time, (c1.chi-w1f.chi)*1e6, 'Color', c(1,:)); hold on
plot(w1s.time, (c1.chi-w1s.chi)*1e6, '--', 'Color', c(1,:))
plot([2200 2200], [0 50], 'k--');
xlim([2000 2300])
ylim([0 50])
ylabel('\Delta CO_2 (ppm)')
title('(b) EW CO_2 concentration reduction')

subplot(3,2,3)
plot(c1.time, c1.DTa, 'k'); hold on
plot(w1f.time, w1f.DTa, 'Color', c(1,:))
plot(w1s.time, w1s.DTa, '--', 'Color', c(1,:))
plot([2200 2200], [0 2], 'k--');
xlim([2000 2300])
ylabel('\Delta SAT (^{\circ}C)')
title('(c) Surface temperature increase')

subplot(3,2,4)
plot(c1.time, c1.DTa-w1f.DTa, 'Color', c(1,:)); hold on
plot(w1s.time, c1.DTa-w1s.DTa, '--', 'Color', c(1,:))
plot([2200 2200], [0 .25], 'k--');
xlim([2000 2300])
ylim([0 0.25])
ylabel('\Delta SAT_{control} - \Delta SAT_{ew} (^{\circ}C)')
title('(d) EW surface temperature reduction')

subplot(3,2,5)
plot(c1.time, c1.pH_carb, 'k'); hold on
plot(w1f.time, w1f.pH_carb, 'Color', c(1,:))
plot(w1s.time, w1s.pH_carb, '--', 'Color', c(1,:))
plot([2200 2200], [7.9 8.1], 'k--');
xlim([2000 2300])
title('(e) Surface average pH')
ylabel('pH')
xlabel('Year')

subplot(3,2,6)
plot(c1.time, -1*(c1.pH_carb-w1f.pH_carb), 'Color', c(1,:)); hold on
plot(w1s.time, -1*(c1.pH_carb-w1s.pH_carb), '--', 'Color', c(1,:))
xlim([2000 2300])
plot([2200 2200], [0 0.06], 'k--');
ylabel('pH_{ew} - pH_{control}')
xlabel('Year')
title('(f) Surface pH increase')
%% Sensitivities of pCO2, T, and pH path dependency among various setups
c2 = load('terr_r26l18t08d05_ctrl.mat');
w2f = load('terr_r26l18t08d05_ewf.mat');
w2s = load('terr_r26l18t08d05_ews.mat');
c3 = load('terr_r45l06t08d05_ctrl.mat');
w3f = load('terr_r45l06t08d05_ewf.mat');
w3s = load('terr_r45l06t08d05_ews.mat');
c4 = load('terr_r45l18t14d08_ctrl.mat');
w4f = load('terr_r45l18t14d08_ewf.mat');
w4s = load('terr_r45l18t14d08_ews.mat');
%%
figure;
subplot(1,3,1)
p(1) = plot(w1s.time, (w1s.chi-w1f.chi)*1e6); hold on
p(2) = plot(w2s.time, (w2s.chi-w2f.chi)*1e6)
p(3) = plot(w3s.time, (w3s.chi-w3f.chi)*1e6)
p(4) = plot(w4s.time, (w4s.chi-w4f.chi)*1e6)
plot([2200 2200], [-5 25], 'k--')
xlim([2100 2300])
title('(a) CO_2 differences')
ylabel('CO_{2, EW-SLOW} - CO_{2, EW-FAST} (ppm)')
legend(p, 'RCP 4.5, high \lambda, weak overturning', ...
    'RCP 2.6, high \lambda, weak overturning', ...
    'RCP 4.5, low \lambda, weak overturning',...
    'RCP 4.5, high \lambda, strong overturning')


subplot(1,3,2)
plot(w1s.time, (c1.DTa-w1s.DTa)./(c1.DTa-w1f.DTa)*100); hold on
plot(w1s.time, (c2.DTa-w2s.DTa)./(c2.DTa-w2f.DTa)*100);
plot(w1s.time, (c3.DTa-w3s.DTa)./(c3.DTa-w3f.DTa)*100);
plot(w1s.time, (c4.DTa-w4s.DTa)./(c4.DTa-w4f.DTa)*100);
plot([2200 2200], [40 110], 'k--')
xlim([2100 2300])
title('(b) Percentage of EW-FAST warming reduction achieved by EW-SLOW')
ylabel('(\Delta SAT_{ctrl} - \Delta SAT_{EW-SLOW})/(\Delta SAT_{ctrl} - \Delta SAT_{EW-FAST}) \times 100')
%legend('89.6', '90.5', '72.9', '91.7')

subplot(1,3,3)
plot(w1s.time, -1*(w1f.pH_carb - w1s.pH_carb)); hold on;
plot(w1s.time, -1*(w2f.pH_carb - w2s.pH_carb));
plot(w1s.time, -1*(w3f.pH_carb - w3s.pH_carb));
plot(w1s.time, -1*(w4f.pH_carb - w4s.pH_carb));
plot([2200 2200], [-0.03 0.01], 'k--')
title('(c) Surface average pH differences')
ylabel('pH_{EW-SLOW} - pH_{EW-FAST}')
xlim([2100, 2300])

%% controls on timing of peak weathering cooling
figure;
plot(c1.time(2:end), diff(c1.N-w1f.N)*365,'k-'); hold on
plot(c1.time(2:end), diff(c1.DR-w1f.DR)*365, 'k--');
plot(c1.time(2:end), diff(c3.N-w3f.N)*365,'b-');
plot(c1.time(2:end), diff(c3.DR-w3f.DR)*365, 'b--')
xlim([2105 2300])
ylabel('W m^{-2} yr^{-1}')
xlabel('Year')
plot([2191 2191], [-3.5e-3, 1.5e-3], 'k:')
legend('d/dt(N_{ctrl}-N_{EW}), standard', ...
    'd/dt(R_{ctrl}-R_{EW}), standard',...
    'd/dt(N_{ctrl}-N_{EW}), low \lambda', ...
    'd/dt(R_{ctrl}-R_{EW}), low \lambda')