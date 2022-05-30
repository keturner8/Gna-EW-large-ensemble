myDir = '/Users/keturner/CE_experiments/lambda_runs/control/';
myFiles = dir(fullfile(myDir,'g2_rcp26_ctrl_*.mat'));

dim = size(myFiles);
dim_files = dim(1);
dim_time = 451;

rcp85_ens = zeros(dim_files, dim_time);

for k = 1:length(myFiles)
    file = load(myFiles(k).name);
    if k==1
        time = file.MM(1:365:end);
    end
    rcp85_ens(k,:) = file.DTa(1:365:end);
    k
end

rcp26_med = median(rcp26_ens, 1);
rcp45_med = median(rcp45_ens, 1);
rcp85_med = median(rcp85_ens, 1);
%%
figure;

subplot(1,3,1)
plot(time, rcp26_ens' - rcp26_ens(:,150)', 'Color', [.5 .5 .5]); hold on
p(1) = plot(time, rcp26_ens(1,:) - rcp26_ens(1,150)', 'Color', [.5 .5 .5])
p(2) = plot(time, rcp26_med - rcp26_med(150), 'k', 'LineWidth', 2)
xlim([2000 2300])
ylim([0 11])
title('RCP 2.6, n=240')
legend(p, 'Individual ensemble members', 'Ensemble median')
ylabel('\Delta T (^{\circ}C)')

subplot(1,3,2)
plot(time, rcp45_ens' - rcp45_ens(:,150)', 'Color', [.5 .5 .5]); hold on
plot(time, rcp45_med - rcp45_med(150), 'k', 'LineWidth', 2)
xlim([2000 2300])
ylim([0 11])
title({'Control surface warming', 'RCP 4.5, n=239'})
xlabel('Year')

subplot(1,3,3)
plot(time, rcp85_ens' - rcp85_ens(:,150)', 'Color', [.5 .5 .5]); hold on
plot(time, rcp85_med - rcp85_med(150), 'k', 'LineWidth', 2)
xlim([2000 2300])
ylim([0 11])
title('RCP 8.5, n=228')