tau = 0.06:0.02:0.14;
del = 0.3:.1:0.8;
rcp_num = [4.5 8.5 2.6]; %2.6 goes at the end because it's the newest set of runs
lambda = 0.6:0.2:2;

empty_mat = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat2 = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat3 = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat4 = nan(length(del), length(tau), length(lambda), length(rcp_num));
empty_mat5 = nan(length(del), length(tau), length(lambda), length(rcp_num));
%% create matrices for peak warming
myFolder = "/Users/keturner/CE_experiments/lambda_runs/control/";

filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
%%
% %% loop for timing of warming
% for i=1:length(theFiles)
%     dummy_ctrl = load(myFolder+theFiles(i).name);
% %    dummy_ccs = load("/Users/keturner/CE_experiments/lambda_runs/CE/CCS/" + ...
% %        theFiles(i).name(1:8)+"_c"+theFiles(i).name(14:end));
%     dummy_ew = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW/" + ...
%         theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end));
%     dummy_ew_slow = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW_half/" + ...
%         theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end-4)+"_half_rate.mat");
%     
%     i/720
%     
% %     if i==1
% %         time_lim = dummy_ctrl.MM>2000;
% %     end
%     
%     idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
%         find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
%     
%     
%     ew_impact_fast = dummy_ctrl.DTa - dummy_ew.DTa;
%     ew_impact_fast_MAX = max(ew_impact_fast);
%     ew_impact_fast_2100 = ew_impact_fast(dummy_ctrl.MM==2100);
%     ew_impact_fast_2200 = ew_impact_fast(dummy_ctrl.MM==2200);
%     
%     ew_impact_slow = dummy_ctrl.DTa - dummy_ew_slow.DTa;
%     ew_impact_slow_MAX = max(ew_impact_slow);
%     %ew_impact_slow_2100 = ew_impact_slow(dummy_ctrl.MM==2100);
%     ew_impact_slow_2200 = ew_impact_slow(dummy_ctrl.MM==2200);
%     
%     empty_mat(idx(1),idx(2),idx(3),idx(4)) = ew_impact_fast_MAX;
%     empty_mat2(idx(1),idx(2),idx(3),idx(4)) = ew_impact_slow_MAX;
%     empty_mat3(idx(1),idx(2),idx(3),idx(4)) = ew_impact_fast_2100;
%     empty_mat4(idx(1),idx(2),idx(3),idx(4)) = ew_impact_fast_2200;
%     empty_mat5(idx(1),idx(2),idx(3),idx(4)) = ew_impact_slow_2200;
%     clear dummy_ctrl dummy_ew dummy_ew_slow %dummy_ccs
% end
%%
%% loop for [H+]
for i=1:length(theFiles)
    dummy_ctrl = load(myFolder+theFiles(i).name);
%    dummy_ccs = load("/Users/keturner/CE_experiments/lambda_runs/CE/CCS/" + ...
%        theFiles(i).name(1:8)+"_c"+theFiles(i).name(14:end));
    dummy_ew = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end));
    dummy_ew_slow = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW_half/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end-4)+"_half_rate.mat");
    
    i/length(theFiles)
    
%     if i==1
%         time_lim = dummy_ctrl.MM>2000;
%     end
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    
    %ew_impact_fast = dummy_ctrl.H_carb - dummy_ew.H_carb;
    %ew_impact_fast_MAX = max(ew_impact_fast);
    %ew_impact_fast_2100 = ew_impact_fast(dummy_ctrl.MM==2100);
    %ew_impact_fast_2200 = ew_impact_fast(dummy_ctrl.MM==2200);
    
    %ew_impact_slow = dummy_ctrl.H_carb - dummy_ew_slow.H_carb;
    %ew_impact_slow_MAX = max(ew_impact_slow);
    %ew_impact_slow_2100 = ew_impact_slow(dummy_ctrl.MM==2100);
    %ew_impact_slow_2200 = ew_impact_slow(dummy_ctrl.MM==2200);
    
    %empty_mat(idx(1),idx(2),idx(3),idx(4)) = ew_impact_fast_MAX;
    %empty_mat2(idx(1),idx(2),idx(3),idx(4)) = ew_impact_slow_MAX;
    empty_mat3(idx(1),idx(2),idx(3),idx(4)) = dummy_ctrl.DTa(dummy_ctrl.MM==2100) / max(dummy_ctrl.DTa);
    empty_mat4(idx(1),idx(2),idx(3),idx(4)) = dummy_ctrl.DTa(dummy_ctrl.MM==2100) - dummy_ew.DTa(dummy_ctrl.MM==2100) ;
    empty_mat2(idx(1),idx(2),idx(3),idx(4)) = max(dummy_ctrl.DTa - dummy_ew.DTa);
    %empty_mat5(idx(1),idx(2),idx(3),idx(4)) = dummy_ew_slow.(dummy_ctrl.MM==2200);
    clear dummy_ctrl dummy_ew dummy_ew_slow %dummy_ccs
end
%%
figure;
for i=1:5
    ZEC_prog = reshape(empty_mat3(:,i,:,3), 48,1);
    CDR_prog = reshape(empty_mat4(:,i,:,3)./empty_mat2(:,i,:,3), 48,1);

    scatter(ZEC_prog, CDR_prog); hold on
end

figure;
scatter(empty_mat3(1,1,1,3), empty_mat4(1,1,1,3)./empty_mat2(1,1,1,3)); hold on
scatter(empty_mat3(1,2,1,3), empty_mat4(1,2,1,3)./empty_mat2(1,2,1,3));
scatter(empty_mat3(1,3,1,3), empty_mat4(1,3,1,3)./empty_mat2(1,3,1,3));
scatter(empty_mat3(1,4,1,3), empty_mat4(1,4,1,3)./empty_mat2(1,4,1,3));
scatter(empty_mat3(1,5,1,3), empty_mat4(1,5,1,3)./empty_mat2(1,5,1,3));
legend

c = zeros(5,3,8);
c(:,:,6) =  [0.9000    0.9000    1.0000;
    0.7000    0.7000    1.0000;
    0.5000    0.5000    1.0000;
    0.3000    0.3000    1.0000;
    0.1000    0.1000    1.0000];
c(:,:,8) = [0.0556    0.0556    0.0816;
    0.2222    0.2222    0.3108;
    0.3889    0.4149    0.5139;
    0.5556    0.6441    0.6806;
    0.7613    0.8472    0.8472];

figure;
for j=1:3
    for i=1:5
     scatter(empty_mat3(3,i,4,j), empty_mat4(3,i,4,j)./empty_mat2(3,i,4,j), [], c(i,:)); hold on 
    end
end
plot([.45 1], [.45 1])

%% loop for overturnings
for i=1:length(theFiles)
    dummy_ctrl = load(myFolder+theFiles(i).name);
%    dummy_ccs = load("/Users/keturner/CE_experiments/lambda_runs/CE/CCS/" + ...
%        theFiles(i).name(1:8)+"_c"+theFiles(i).name(14:end));
    
    i/length(theFiles)
    
%     if i==1
%         time_lim = dummy_ctrl.MM>2000;
%     end
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    
    circ_pi = dummy_ctrl.VTna(dummy_ctrl.MM==1860);
    empty_mat(idx(1),idx(2),idx(3),idx(4)) = circ_pi;

    clear dummy_ctrl
end

%% loop for CO2 differences
for i=1:length(theFiles)
    dummy_ctrl = load(myFolder+theFiles(i).name);
    dummy_ew = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end));
    dummy_ew_slow = load("/Users/keturner/CE_experiments/lambda_runs/CE/EW_half/" + ...
        theFiles(i).name(1:8)+"_2w2000"+theFiles(i).name(14:end-4)+"_half_rate.mat");
    
    i/length(theFiles)
    
    
    idx = [find(round(del,1)==round(dummy_ctrl.FRAC,1)), find(tau==dummy_ctrl.tau), ...
        find(round(lambda,1) == round(dummy_ctrl.new_lambda,1)), find(rcp_num==str2num(theFiles(i).name(7:8))/10)];
    
    
    co2_ctrl = dummy_ctrl.chi(dummy_ctrl.MM==2200);
    co2_ew = dummy_ew.chi(dummy_ctrl.MM==2200);
    co2_ew_slow = dummy_ew_slow.chi(dummy_ctrl.MM==2200);
    empty_mat(idx(1),idx(2),idx(3),idx(4)) = co2_ctrl;
    empty_mat2(idx(1),idx(2),idx(3),idx(4)) = co2_ew;
    empty_mat3(idx(1),idx(2),idx(3),idx(4)) = co2_ew_slow;

    clear dummy_ctrl dummy_ew dummy_ew_slow
end

