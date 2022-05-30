
directory = "/Users/keturner/CE_experiments/lambda_runs/";

ex1_ctrl = load(directory+"control/g2_rcp26_ctrl_d60_t080_l18.mat");
ex1_fullw = load(directory+"CE/EW/g2_rcp26_2w2000_d60_t080_l18.mat");
ex1_fullc = load(directory+"CE/CCS/g2_rcp26_c_d60_t080_l18.mat");
ex1_halfw = load(directory+"CE/EW_half/g2_rcp26_2w2000_d60_t080_l18_half_rate.mat");
ex1_halfc = load(directory+"CE/CCS_half/g2_rcp26_2c2000_d60_t080_l18_half_rate.mat");


ex2_ctrl = load(directory+"control/g2_rcp85_ctrl_d30_t120_l06.mat");
ex2_fullw = load(directory+"CE/EW/g2_rcp85_2w2000_d30_t120_l06.mat");
ex2_fullc = load(directory+"CE/CCS/g2_rcp85_c_d30_t120_l06.mat");
ex2_halfw = load(directory+"CE/EW_half/g2_rcp85_2w2000_d30_t120_l06_half_rate.mat");
ex2_halfc = load(directory+"CE/CCS_half/g2_rcp85_2c2000_d30_t120_l06_half_rate.mat");

tt = ex1_ctrl.MM;
%% Figure 1
c = get(gca, 'colororder');

figure;
subplot(2,2,1)
p(1) = plot(tt, ex1_ctrl.chi * 1e6, 'k'); hold on
p(2) = plot(tt, ex2_ctrl.chi * 1e6, 'Color', [.5 .5 .5]);
plot(tt, ex1_fullw.chi * 1e6, 'Color', c(1,:));
plot(tt, ex2_fullw.chi * 1e6, 'Color', c(2,:));
plot(tt, ex1_halfw.chi * 1e6, '--', 'Color', c(1,:));
plot(tt, ex2_halfw.chi * 1e6, '--', 'Color', c(2,:));
plot([2200 2200], [400 1400], 'k--');
xlim([2000 2300])
title('(a) Atmospheric CO_2 concentrations')
ylabel('CO_2 (ppm)')
xlabel('Year')
legend(p, 'Control, RCP 2.6, high \lambda', 'Control, RCP 8.5, low \lambda')

subplot(2,2,2)
plot(tt, (ex1_ctrl.chi - ex1_fullw.chi) * 1e6, 'Color', c(1,:)); hold on
plot(tt, (ex1_ctrl.chi - ex1_halfw.chi) * 1e6, '--', 'Color', c(1,:));
plot(tt, (ex2_ctrl.chi - ex2_fullw.chi) * 1e6, 'Color', c(2,:));;
plot(tt, (ex2_ctrl.chi - ex2_halfw.chi) * 1e6, '--', 'Color', c(2,:));
plot([2200 2200], [0 80], 'k--');
xlim([2000 2300])
ylim([0 80])
title('(b) EW CO_2 concentration reduction')
ylabel('\Delta CO_2 (ppm)')
xlabel('Year')
legend('2 PgC/yr EW \times 100 yr, RCP 2.6 & high \lambda', ...
    '1 PgC/yr EW \times 200 yr, RCP 2.6 & high \lambda', ...
    '2 PgC/yr EW \times 100 yr, RCP 8.5 & low \lambda',...
    '1 PgC/yr EW \times 200 yr, RCP 8.5 & low \lambda')

subplot(2,2,3)
plot(tt, ex1_ctrl.DTa, 'k'); hold on
plot(tt, ex2_ctrl.DTa, 'k');
plot(tt, ex1_fullw.DTa, 'Color', c(1,:));
plot(tt, ex2_fullw.DTa, 'Color', c(2,:));
plot(tt, ex1_halfw.DTa, '--', 'Color', c(1,:));
plot(tt, ex2_halfw.DTa, '--', 'Color', c(2,:));
plot([2200 2200], [0 12], 'k--');
xlim([2000 2300])
title('(c) Surface temperature increase')
ylabel('\Delta SAT (^{\circ}C)')
xlabel('Year')

subplot(2,2,4)
plot(tt, ex1_ctrl.DTa - ex1_fullw.DTa, 'Color', c(1,:)); hold on
plot(tt, ex1_ctrl.DTa - ex1_halfw.DTa, '--', 'Color', c(1,:));
plot(tt, ex2_ctrl.DTa - ex2_fullw.DTa, 'Color', c(2,:));
plot(tt, ex2_ctrl.DTa - ex2_halfw.DTa, '--', 'Color', c(2,:));
xlim([2000 2300])
plot([2200 2200], [0 0.6], 'k--');
%ylim([0 80])
title('(d) EW surface temperature reduction')
ylabel('\Delta SAT_{control} - \Delta SAT_{ew} (^{\circ}C)')
xlabel('Year')

%% pH
c = get(gca, 'colororder');
subplot(2,2,1)
q(1) = plot(tt, ex1_ctrl.Iocean - ex1_ctrl.Iocean(1), 'k'); hold on
q(2) = plot(tt, ex2_ctrl.Iocean - ex1_ctrl.Iocean(1),  'Color', [.5 .5 .5]);
q(3) = plot(tt, ex1_fullw.Iocean - ex1_ctrl.Iocean(1), 'Color', c(1,:));
q(4) = plot(tt, ex2_fullw.Iocean - ex1_ctrl.Iocean(1), 'Color', c(2,:));
plot(tt, ex1_fullc.Iocean - ex1_ctrl.Iocean(1), '-.', 'Color', c(1,:));
plot(tt, ex2_fullc.Iocean - ex1_ctrl.Iocean(1), '-.', 'Color', c(2,:));
plot(tt, ex1_halfw.Iocean - ex1_ctrl.Iocean(1), '--', 'Color', c(1,:));
plot(tt, ex2_halfw.Iocean - ex1_ctrl.Iocean(1), '--', 'Color', c(2,:));
plot(tt, ex1_halfc.Iocean - ex1_ctrl.Iocean(1), ':', 'Color', c(1,:));
plot(tt, ex2_halfc.Iocean - ex1_ctrl.Iocean(1), ':', 'Color', c(2,:));
plot([2200 2200], [100 800], 'k--', 'LineWidth', 1)
legend(q, 'Control, RCP 2.6, high \lambda', 'Control, RCP 8.5, low \lambda',...
    'CDR, RCP 2.6, high \lambda', 'CDR, RCP 8.5, low \lambda')
xlim([2000 2300])
title('(a) Changes in ocean carbon inventory')
ylabel('\Delta I_{ocean} (PgC)')

subplot(2,2,2)
plot(tt, ex1_fullw.Iocean - ex1_ctrl.Iocean, 'Color', c(1,:)); hold on
plot(tt, ex2_fullw.Iocean - ex2_ctrl.Iocean, 'Color', c(2,:));
plot(tt, ex1_halfw.Iocean - ex1_ctrl.Iocean, '--', 'Color', c(1,:));
plot(tt, ex2_halfw.Iocean - ex2_ctrl.Iocean, '--', 'Color', c(2,:));
plot(tt, ex1_fullc.Iocean - ex1_ctrl.Iocean, '-.', 'Color', c(1,:));
plot(tt, ex2_fullc.Iocean - ex2_ctrl.Iocean, '-.', 'Color', c(2,:));
plot(tt, ex1_halfc.Iocean - ex1_ctrl.Iocean, ':', 'Color', c(1,:));
plot(tt, ex2_halfc.Iocean - ex2_ctrl.Iocean, ':', 'Color', c(2,:));
plot([2200 2200], [-100 180], 'k--', 'LineWidth', 1)
ylim([-100 180])
xlim([2000 2300])
title('(c) Ocean carbon changes from CDR')
ylabel('I_{ocean, CDR} - I_{ocean, control} (PgC)')

subplot(2,2,3)
plot(tt, ex1_ctrl.pH_carb, 'k'); hold on
plot(tt, ex2_ctrl.pH_carb, 'Color', [.5 .5 .5]);
plot(tt, ex1_fullw.pH_carb, 'Color', c(1,:));
plot(tt, ex2_fullw.pH_carb, 'Color', c(2,:));
plot(tt, ex1_fullc.pH_carb, '-.', 'Color', c(1,:));
plot(tt, ex2_fullc.pH_carb, '-.', 'Color', c(2,:));
plot(tt, ex1_halfw.pH_carb, '--', 'Color', c(1,:));
plot(tt, ex2_halfw.pH_carb, '--', 'Color', c(2,:));
plot(tt, ex1_halfc.pH_carb, ':', 'Color', c(1,:));
plot(tt, ex2_halfc.pH_carb, ':', 'Color', c(2,:));
plot([2200 2200], [7.6 8.1], 'k--');
xlim([2000 2300])
title('(c) Surface average pH')
ylabel('pH')
xlabel('Year')


subplot(2,2,4)
plot(tt, ex1_fullw.pH_carb - ex1_ctrl.pH_carb, 'Color', c(1,:)); hold on
plot(tt, ex1_halfw.pH_carb - ex1_ctrl.pH_carb, '--', 'Color', c(1,:));
plot(tt, ex1_fullc.pH_carb - ex1_ctrl.pH_carb, '-.', 'Color', c(1,:));
plot(tt, ex1_halfc.pH_carb - ex1_ctrl.pH_carb, ':', 'Color', c(1,:));
plot(tt, ex2_fullw.pH_carb - ex2_ctrl.pH_carb, 'Color', c(2,:));
plot(tt, ex2_halfw.pH_carb - ex2_ctrl.pH_carb, '--', 'Color', c(2,:));
plot(tt, ex2_fullc.pH_carb - ex2_ctrl.pH_carb, '-.', 'Color', c(2,:));
plot(tt, ex2_halfc.pH_carb - ex2_ctrl.pH_carb, ':', 'Color', c(2,:));
xlim([2000 2300])
plot([2200 2200], [-0.01 0.08], 'k--');
%ylim([0 80])
title('(d) Surface pH increase')
ylabel('pH_{EWl} - pH_{control}')
xlabel('Year')
legend('2 PgC/yr EW \times 100 yr, RCP 2.6 & high \lambda', ...
    '1 PgC/yr EW \times 200 yr, RCP 2.6 & high \lambda', ...
    '2 PgC/yr EW \times 100 yr, RCP 8.5 & low \lambda',...
    '1 PgC/yr EW \times 200 yr, RCP 8.5 & low \lambda')