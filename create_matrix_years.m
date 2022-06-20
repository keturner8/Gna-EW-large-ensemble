tau = 0.06:0.02:0.14;
del = 0.3:.1:0.8;
rcp_num = [4.5 8.5 2.6]; %2.6 goes at the end because it's the newest set of runs
lambda = 0.6:0.2:2.0;

empty_mat = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat2 = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat3 = nan(length(del), length(tau), length(lambda), length(rcp_num));
%% create matrices for peak warming
myFolder = "/Users/keturner/CE_experiments/lambda_runs/control/";

filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

%% loop for timing of warming
for i=1:length(theFiles)
    dummy_ctrl = load(myFolder+theFiles(i).name);
    dummy_ccs = load("/Users/keturner/CE_experiments/lambda_runs/CE/CCS/" + ...
        theFiles(i).name(1:8)+"_c"+theFiles(i).name(14:end));
    dummy_ew = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end));
    
    i/720
    
    if i==1
        time_lim = dummy_ctrl.MM>2000;
    end
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    warming_peak_times_ctrl = islocalmax(dummy_ctrl.DTa) & time_lim;
    warming_peak_times_ccs = islocalmax(dummy_ctrl.DTa - dummy_ccs.DTa) & time_lim;
    warming_peak_times_ew = islocalmax(dummy_ctrl.DTa - dummy_ew.DTa) & time_lim;
    
    if sum(warming_peak_times_ctrl) > 1
        times = find(warming_peak_times_ctrl);
        peak_yr_ctrl = dummy_ctrl.MM(times(1));
    elseif sum(warming_peak_times_ctrl) == 0
        peak_yr_ctrl = 2301;
    else
        peak_yr_ctrl = dummy_ctrl.MM(warming_peak_times_ctrl);
    end
    
    if sum(warming_peak_times_ccs) > 1
        times = find(warming_peak_times_ccs);
        peak_yr_ccs = dummy_ctrl.MM(times(1));% pick one of these times to be representative, say the first one?
    elseif sum(warming_peak_times_ccs) == 0
        peak_yr_ccs = 2301; % peak time is after 2300, so would need a longer run, which we can do later...
    else
        peak_yr_ccs = dummy_ctrl.MM(warming_peak_times_ccs);
        %peak_temp = dummy_ctrl.DTa(warming_peak_times);% - dummy_ce.DTa(warming_peak_times);
    end
    
    if sum(warming_peak_times_ew) > 1
        times = find(warming_peak_times_ew);
        peak_yr_ew = dummy_ctrl.MM(times(1));% pick one of these times to be representative, say the first one?
    elseif sum(warming_peak_times_ew) == 0
        peak_yr_ew = 2301; % peak time is after 2300, so would need a longer run, which we can do later...
    else
        peak_yr_ew = dummy_ctrl.MM(warming_peak_times_ew);
        %peak_temp = dummy_ctrl.DTa(warming_peak_times);% - dummy_ce.DTa(warming_peak_times);
    end
    
    empty_mat(idx(1),idx(2),idx(3),idx(4)) = peak_yr_ctrl;
    empty_mat2(idx(1),idx(2),idx(3),idx(4)) = peak_yr_ccs;
    empty_mat3(idx(1),idx(2),idx(3),idx(4)) = peak_yr_ew;
    clear dummy_ctrl dummy_ccs dummy_ew
end

%% loop for warming at 2100 when emissions end

for i=1:length(theFiles)
    dummy_ctrl = load(theFiles(i).name);
    
    i
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(lambda == dummy_ctrl.new_lambda), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    if i==1 % only have to do this once
        idx_2100 = dummy_ctrl.MM == 2100;
    end
    
    empty_mat3(idx(1),idx(2),idx(3),idx(4)) = dummy_ctrl.DTa(idx_2100);
    clear dummy
end

%% loop for calculating when EW methods are approximately equivalent

for i=1:length(theFiles)
    dummy_ctrl = load(myFolder+theFiles(i).name);
    dummy_ew1 = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end));
    dummy_ew2 = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW_half/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end-4) + "_half_rate.mat");
    
    i/720
    
    if i==1
        time_lim = dummy_ctrl.MM>2100;
        time_subsamp = dummy_ctrl.MM(time_lim);
        eps = 0.95; %setting percentage to which slow EW is within fast EW impact
    end
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    impact_ew1 = dummy_ctrl.DTa(time_lim) - dummy_ew1.DTa(time_lim);
    impact_ew2 = dummy_ctrl.DTa(time_lim) - dummy_ew2.DTa(time_lim);
    
    ratio = impact_ew2./impact_ew1;
    equiv = ratio >= eps;
    
    if sum(equiv) == 0
        a = nan;
    else
        a = min(time_subsamp(equiv));
    end
    
    empty_mat(idx(1),idx(2),idx(3),idx(4)) = a;
    clear dummy_ctrl dummy_ew1 dummy_ew2
end