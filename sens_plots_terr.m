c = dir('terr*ctrl.mat');
w1 = dir('terr*ewf.mat');
w2 = dir('terr*ews.mat');

lambda = 0.6:0.2:2.0;
cb = cmocean('matter');

load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/pco2_atm_2200.mat');
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/adjustment_slices.mat');
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/pH_2200.mat');
figure;
for k=1:length(lambda)
    A1 = (co2_ctrl(:,:,k,3) - co2_ews(:,:,k,3)) ./ (co2_ctrl(:,:,k,3) - co2_ewf(:,:,k,3));
    A2 = (co2_ctrl(:,:,k,1) - co2_ews(:,:,k,1)) ./ (co2_ctrl(:,:,k,1) - co2_ewf(:,:,k,1));

    subplot(1,3,1)
    scatter(reshape(A1, numel(A1), 1), repmat(lambda(k), 1, numel(A1)), [], [.6 .6 .6], 'd'); hold on
    scatter(reshape(A2, numel(A2), 1), repmat(lambda(k), 1, numel(A2)), [], [.6 .6 .6]);

    B1 = ew_impact_slow_2200(:,:,k,3) ./ ew_impact_fast_2200(:,:,k,3);
    B2 = ew_impact_slow_2200(:,:,k,1) ./ ew_impact_fast_2200(:,:,k,1);

    subplot(1,3,2)
    scatter(reshape(B1, numel(B1), 1), repmat(lambda(k), 1, numel(B1)), [], [.6 .6 .6], 'd'); hold on
    scatter(reshape(B2, numel(B2), 1), repmat(lambda(k), 1, numel(B2)), [], [.6 .6 .6]);

    C1 = pH_ewf_2200(:,:,k,3) - pH_ews_2200(:,:,k,3);
    C2 = pH_ewf_2200(:,:,k,1) - pH_ews_2200(:,:,k,1);

    subplot(1,3,3)
    scatter(reshape(C1, numel(C1), 1), repmat(lambda(k), 1, numel(C1)), [], [.6 .6 .6], 'd'); hold on
    scatter(reshape(C2, numel(C2), 1), repmat(lambda(k), 1, numel(C2)), [], [.6 .6 .6]);
end

for i = 1:11
    f_ctrl = load(c(i).name);
    f_ewf = load(w1(i).name);
    f_ews = load(w2(i).name);

    if i==1
        idx_2200 = find(f_ctrl.time == 2200);
    end

    if sum(c(i).name(7:8) == '45')
        mkr_sty = 'o';
    else
        mkr_sty = 'd';
    end

    if sum(c(i).name(13:14) == '14')
        mkr_clr = cb(200,:);

    else
        mkr_clr = cb(20,:);
    end

    subplot(1,3,1)
    scatter((f_ctrl.chi(idx_2200) - f_ews.chi(idx_2200)) / (f_ctrl.chi(idx_2200) - f_ewf.chi(idx_2200)),...
        f_ewf.lambda, [], mkr_clr, mkr_sty, 'filled');
    hold on;

    subplot(1,3,2)
    scatter((f_ctrl.DTa(idx_2200) - f_ews.DTa(idx_2200)) / (f_ctrl.DTa(idx_2200) - f_ewf.DTa(idx_2200)), ...
        f_ewf.lambda, [], mkr_clr, mkr_sty, 'filled');
    hold on;

    subplot(1,3,3)
    scatter(f_ewf.pH_carb(idx_2200) - f_ews.pH_carb(idx_2200), f_ewf.lambda, ...
        [], mkr_clr, mkr_sty, 'filled');
    hold on;
end

subplot(1,3,1)
f=get(gca,'Children');
legend([f(7), f(8), f(2), f(1), f(end), f(end-1)],'RCP 2.6, weak overturning', ...
    'RCP 2.6, strong overturning', 'RCP 4.5, weak overturning', ...
    'RCP 4.5, strong overturning', 'RCP 2.6, no land setup', ...
    'RCP 4.5, no land setup')
xlim([.9 1.1])
title('(a) Relative atmospheric CO_2 difference');
ylabel('Climate feedback parameter \lambda, W m^{-2} K^{-1}');
xlabel('(pCO_{2,control} - pCO_{2,EW-SLOW})/(pCO_{2,control} - pCO_{2,EW-FAST}) |_{t=2200}')
box on
ylim([.5 2.1])
subplot(1,3,2)
title('(b) Relative surface temperature difference');
box on
ylim([.5 2.1])
xlabel('(\Delta T_{control} - \Delta T_{EW-SLOW})/(\Delta T_{control} - \Delta T_{EW-FAST})|_{t=2200}')
subplot(1,3,3)
title('(c) Surface ocean pH difference');
box on
ylim([.5 2.1])
xlabel('(pH_{EW-FAST} - pH_{EW_SLOW})|_{t=2200}')