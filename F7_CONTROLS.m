directory = "/Users/keturner/CE_experiments/lambda_runs/" ;

tau = [0.06 0.14];
del = [0.3 0.8];
rcp_num = [8.5 2.6];
lambda = [0.6 1.8];

std_run_ctrl = load(directory+"control/g2_rcp45_ctrl_d50_t080_l18.mat");
std_run_weat = load(directory+"CE/EW/g2_rcp45_2w2000_d50_t080_l18.mat");

llow_run_ctrl = load(directory+"control/g2_rcp45_ctrl_d50_t080_l06.mat");
llow_run_weat = load(directory+"CE/EW/g2_rcp45_2w2000_d50_t080_l06.mat");

Rhigh_run_ctrl = load(directory+"control/g2_rcp26_ctrl_d50_t080_l18.mat");
Rhigh_run_weat = load(directory+"CE/EW/g2_rcp26_2w2000_d50_t080_l18.mat");

thigh_run_ctrl = load(directory+"control/g2_rcp45_ctrl_d50_t140_l18.mat");
thigh_run_weat = load(directory+"CE/EW/g2_rcp45_2w2000_d50_t140_l18.mat");

dhigh_run_ctrl = load(directory+"control/g2_rcp45_ctrl_d80_t080_l18.mat");
dhigh_run_weat = load(directory+"CE/EW/g2_rcp45_2w2000_d80_t080_l18.mat");

tt = std_run_ctrl.MM;
tt_subsamp = tt>=2100;
tt_subsamp2 = tt(tt_subsamp);
tt_subsamp3 = tt_subsamp2(2:end);
%%
test_std = abs(diff(std_run_ctrl.N(tt_subsamp) - std_run_weat.N(tt_subsamp)) - ...
    diff(std_run_ctrl.DR(tt_subsamp) - std_run_weat.DR(tt_subsamp)));

yr_std = tt_subsamp3(test_std == min(test_std));
%
test_llow = abs(diff(llow_run_ctrl.N(tt_subsamp) - llow_run_weat.N(tt_subsamp)) - ...
    diff(llow_run_ctrl.DR(tt_subsamp) - llow_run_weat.DR(tt_subsamp)));

yr_llow = tt_subsamp3(test_llow == min(test_llow));
%
test_Rhigh = abs(diff(Rhigh_run_ctrl.N(tt_subsamp) - Rhigh_run_weat.N(tt_subsamp)) - ...
    diff(Rhigh_run_ctrl.DR(tt_subsamp) - Rhigh_run_weat.DR(tt_subsamp)));

yr_Rhigh = tt_subsamp3(test_Rhigh == min(test_Rhigh));
%
test_thigh = abs(diff(thigh_run_ctrl.N(tt_subsamp) - thigh_run_weat.N(tt_subsamp)) - ...
    diff(thigh_run_ctrl.DR(tt_subsamp) - thigh_run_weat.DR(tt_subsamp)));

yr_thigh = tt_subsamp3(test_thigh == min(test_thigh));
%
test_dhigh = abs(diff(dhigh_run_ctrl.N(tt_subsamp) - dhigh_run_weat.N(tt_subsamp)) - ...
    diff(dhigh_run_ctrl.DR(tt_subsamp) - dhigh_run_weat.DR(tt_subsamp)));

yr_dhigh = tt_subsamp3(test_dhigh == min(test_dhigh));
%%
figure;
subplot(3,2,1)
plot(tt, std_run_ctrl.N, 'k'); hold on
plot(tt, std_run_weat.N, 'k:');
xlim([2100 2300])
ylabel('N (W m^{-2})')
title('(a) Heat uptake post-emissions')
legend('Sample control run', 'Sample weathering run')

subplot(3,2,2)
plot(tt, std_run_ctrl.DR, 'k--'); hold on
plot(tt, std_run_weat.DR, 'k-.');
xlim([2100 2300])
title('(b) Radiative forcing post-emissions')
ylabel('\Delta R (W m^{-2})')
legend('Sample control run', 'Sample weathering run')

subplot(3,2,3)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k'); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--')
plot(tt(2:end),365*diff(llow_run_ctrl.N - llow_run_weat.N));
plot(tt(2:end),365*diff(llow_run_ctrl.DR - llow_run_weat.DR))
xlim([2100 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:')
plot([yr_llow yr_llow], [-2e-3, 1e-3], 'b:')
ylim([-2e-3, 1e-3])
ylabel('W m^{-2} yr^{-1}')
title('(c) Varying \lambda')
legend('d/dt(N_{ctrl}-N_{EW}), standard', ...
    'd/dt(R_{ctrl}-R_{EW}), standard',...
    'd/dt(N_{ctrl}-N_{EW}), low \lambda',...
    'd/dt(R_{ctrl}-R_{EW}), low \lambda')

subplot(3,2,4)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k'); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--')
p(1) = plot(tt(2:end),365*diff(Rhigh_run_ctrl.N - Rhigh_run_weat.N));
p(2) = plot(tt(2:end),365*diff(Rhigh_run_ctrl.DR - Rhigh_run_weat.DR))
xlim([2100 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:')
ylabel('W m^{-2} yr^{-1}')
%plot([yr_Rhigh yr_Rhigh], [-2e-3, 1e-3], 'b:')
ylim([-2e-3, 1e-3])
title('(d) Varying emissions')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), RCP 8.5',...
    'd/dt(R_{ctrl}-R_{EW}), RCP 8.5')

subplot(3,2,5)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k'); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--')
p(1)=plot(tt(2:end),365*diff(thigh_run_ctrl.N - thigh_run_weat.N));
p(2)=plot(tt(2:end),365*diff(thigh_run_ctrl.DR - thigh_run_weat.DR))
xlim([2100 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:')
plot([yr_thigh yr_thigh], [-2e-3, 1e-3], 'b:')
ylabel('W m^{-2} yr^{-1}')
ylim([-2e-3, 1e-3])
title('(e) Varying \tau')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), high \tau',...
    'd/dt(R_{ctrl}-R_{EW}), high \tau')
xlabel('Year')

subplot(3,2,6)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k'); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--')
p(1)=plot(tt(2:end),365*diff(dhigh_run_ctrl.N - dhigh_run_weat.N));
p(2)=plot(tt(2:end),365*diff(dhigh_run_ctrl.DR - dhigh_run_weat.DR))
xlim([2100 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:')
plot([yr_dhigh yr_dhigh], [-2e-3, 1e-3], 'b:')
ylabel('W m^{-2} yr^{-1}')
ylim([-2e-3, 1e-3])
title('(f) Varying \delta')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), high \delta',...
    'd/dt(R_{ctrl}-R_{EW}), high \delta')
xlabel('Year')

%%
C = colororder;

subplot(2,2,1)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k', 'LineWidth', 1.5); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--', 'LineWidth', 1.5)
plot(tt(2:end),365*diff(llow_run_ctrl.N - llow_run_weat.N), 'Color',C(1,:), 'LineWidth', 1.5);
plot(tt(2:end),365*diff(llow_run_ctrl.DR - llow_run_weat.DR), '--', 'Color', C(1,:), 'LineWidth', 1.5)
xlim([2101 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:', 'LineWidth', 1.5)
plot([yr_llow yr_llow], [-2e-3, 1e-3], ':', 'Color', C(1,:), 'LineWidth', 1.5)
ylim([-2e-3, 1e-3])
ylabel('W m^{-2} yr^{-1}')
title('(a) Varying \lambda')
legend('d/dt(N_{ctrl}-N_{EW}), standard', ...
    'd/dt(R_{ctrl}-R_{EW}), standard',...
    'd/dt(N_{ctrl}-N_{EW}), low \lambda',...
    'd/dt(R_{ctrl}-R_{EW}), low \lambda')
xlabel('year')

subplot(2,2,4)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k', 'LineWidth', 1.5); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--', 'LineWidth', 1.5)
p(1) = plot(tt(2:end),365*diff(Rhigh_run_ctrl.N - Rhigh_run_weat.N), 'Color', C(2,:), 'LineWidth', 1.5);
p(2) = plot(tt(2:end),365*diff(Rhigh_run_ctrl.DR - Rhigh_run_weat.DR), '--', 'Color', C(2,:), 'LineWidth', 1.5)
xlim([2101 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:', 'LineWidth', 1.5)
ylabel('W m^{-2} yr^{-1}')
plot([yr_Rhigh yr_Rhigh], [-2e-3, 1e-3], ':', 'Color', C(2,:), 'LineWidth', 1.5)
ylim([-2e-3, 1e-3])
title('(d) Varying emissions')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), RCP 2.6',...
    'd/dt(R_{ctrl}-R_{EW}), RCP 2.6')
xlabel('Year')

subplot(2,2,2)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k', 'LineWidth', 1.5); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--', 'LineWidth', 1.5)
p(1)=plot(tt(2:end),365*diff(thigh_run_ctrl.N - thigh_run_weat.N), 'Color', C(3,:), 'LineWidth', 1.5);
p(2)=plot(tt(2:end),365*diff(thigh_run_ctrl.DR - thigh_run_weat.DR), '--', 'Color', C(3,:), 'LineWidth', 1.5)
xlim([2101 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:', 'LineWidth', 1.5)
plot([yr_thigh yr_thigh], [-2e-3, 1e-3], ':' , 'Color', C(3,:), 'LineWidth', 1.5)
ylabel('W m^{-2} yr^{-1}')
ylim([-2e-3, 1e-3])
title('(b) Varying \tau')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), high \tau',...
    'd/dt(R_{ctrl}-R_{EW}), high \tau')
xlabel('Year')

subplot(2,2,3)
plot(tt(2:end),365*diff(std_run_ctrl.N - std_run_weat.N), 'k', 'LineWidth', 1.5); hold on;
plot(tt(2:end),365*diff(std_run_ctrl.DR - std_run_weat.DR), 'k--', 'LineWidth', 1.5)
p(1)=plot(tt(2:end),365*diff(dhigh_run_ctrl.N - dhigh_run_weat.N), 'Color', C(4,:), 'LineWidth', 1.5);
p(2)=plot(tt(2:end),365*diff(dhigh_run_ctrl.DR - dhigh_run_weat.DR), '--', 'Color', C(4,:), 'LineWidth', 1.5)
xlim([2101 2300])
plot([yr_std yr_std], [-2e-3, 1e-3], 'k:', 'LineWidth', 1.5)
plot([yr_dhigh yr_dhigh], [-2e-3, 1e-3], ':', 'Color', C(4,:), 'LineWidth', 1.5)
ylabel('W m^{-2} yr^{-1}')
ylim([-2e-3, 1e-3])
title('(c) Varying \delta')
legend(p, 'd/dt(N_{ctrl}-N_{EW}), high \delta',...
    'd/dt(R_{ctrl}-R_{EW}), high \delta')
xlabel('Year')