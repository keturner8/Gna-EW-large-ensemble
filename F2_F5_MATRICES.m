tau = 0.06:0.02:0.14;
del = 0.3:.1:0.8;
rcp_num = [4.5 8.5 2.6];
lambda = 0.6:0.2:2.0;

% climate feedback parameters for cmip5 vs cmip6 models
lambda_cmip5 = -1*[0.76 0.81 0.92 1.18 1.14 0.64 1.03 0.84 0.91 0.76 1.23 ...
    1.37 1.66 1.76 0.63 0.75 0.80 1.02 0.92 1.53 1.13 1.19 1.23 1.22 1.11 ...
    1.15 1.19 1.43];

lambda_cmip6 = -1*[1.03 0.92 1.82 0.63 0.71 0.74 0.92 0.62 0.65 0.63 0.81 ...
    0.78 1.40 0.82 1.45 1.13 0.63 1.48 0.75 1.54 1.40 1.22 1.10 0.76 1.34 ...
    1.05 0.67];

%% lambda snapshots
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/adjustment_slices.mat');
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/circ_VTna_pi.mat');
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/pH_2200.mat');
load('/Users/keturner/CE_experiments/lambda_runs/matrices_inc2.6/peak_yrs_2.mat');
circ_VTna_pi = squeeze(circ_VTna_pi(:,:,:,1));
%%
figure;
c = get(gca, 'ColorOrder')
for k=1:length(lambda)
    A = ew_impact_fast_MAX(:,:,k,3);
    B1 = ew_impact_fast_2100(:,:,k,3);
    B2 = ew_impact_fast_2200(:,:,k,3);
    C = ew_impact_slow_MAX(:,:,k,3);
    D = ew_impact_slow_2200(:,:,k,3);

    subplot(2,3,1)
    q1 = scatter(100*reshape(B1./A, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    xlim([40 90])
    ylim([.5 2.1])
    if k==1
        title('(a) RCP 2.6')
        ylabel({'Climate feedback parameter \lambda,', 'W m^{-2} K^{-1}'})
    end
    
    subplot(2,3,4)
    q1 = scatter(100*reshape(D./B2, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    xlim([70 95])
    ylim([.5 2.1])
    if k==1
        title('(d) RCP 2.6')
        ylabel({'Climate feedback parameter \lambda,', 'W m^{-2} K^{-1}'})
    end

  
    A = ew_impact_fast_MAX(:,:,k,1);
    B1 = ew_impact_fast_2100(:,:,k,1);
    B2 = ew_impact_fast_2200(:,:,k,1);
    C = ew_impact_slow_MAX(:,:,k,1);
    D = ew_impact_slow_2200(:,:,k,1);

    subplot(2,3,2)
    q1 = scatter(100*reshape(B1./A, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    xlim([40 90])
    ylim([.5 2.1])
    if k==1
        xlabel({'\bf Fraction of peak EW-FAST mitigation achieved by 2100',...\
            '\rm \Delta T_{EW-FAST}(t=2100)/max(\Delta T_{EW-FAST})'})
        title('(b) RCP 4.5')
    end
    
    subplot(2,3,5)
    scatter(100*reshape(D./B2, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    ylim([.5 2.1])
    xlim([70 95])
    if k==1
        xlabel({'\bf Fraction of EW-FAST mitigation achieved by EW-SLOW at 2200',...\
            '\rm \Delta T_{EW-SLOW}(t=2200)/\Delta T_{EW-FAST}(t=2200)'})
        title('(e) RCP 4.5')
    end
 
    A = ew_impact_fast_MAX(:,:,k,2);
    B1 = ew_impact_fast_2100(:,:,k,2);
    B2 = ew_impact_fast_2200(:,:,k,2);
    C = ew_impact_slow_MAX(:,:,k,2);
    D = ew_impact_slow_2200(:,:,k,2);

    idx_locmax = double(peak_yr_ew_2(:,:,k,2)~=2301);
    idx_locmax(idx_locmax==0)=nan;

    subplot(2,3,3)
    q1 = scatter(100*reshape(B1./A.*idx_locmax, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    xlim([40 90])
    ylim([.5 2.1])
    if k==1
        title('(c) RCP 8.5')
    end
    
    subplot(2,3,6)
    scatter(100*reshape(D./B2, numel(A), 1), repmat(lambda(k), 1, numel(A)), [], reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1)); hold on
    cmocean('matter')
    ylim([.5 2.1])
    xlim([70 95])
    if k==1
        title('(f) RCP 8.5')
    end
end
    %legend(p, '\Delta SAT, EW @ 2 PgC/yr x 100 yr', ...
    %    '\Delta SAT, EW @ 1 PgC/yr x 200 yr');
%%
figure;
subplot(1,2,1)
hist(-1*lambda_cmip5, [0.6:0.2:1.8]);
view(90,-90)
xlim([0.5, 1.9])
ylim([0 11])
title({'CMIP5 \lambda,', 'N = 28'})
xlabel('\lambda, W m^{-2} K^{-1}')
ylabel('Frequency (#)')

subplot(1,2,2)
hist(-1*lambda_cmip6, [0.6:0.2:2]);
view(90,-90)
xlim([0.5, 1.9])
ylim([0 11])
title({'CMIP6 \lambda,', 'N = 27'})
ylabel('Frequency (#)')

%% pH figures - differences
lambda_pH = 0.6:0.2:2.0;

figure;
for k=1:length(lambda)
    subplot(1,3,1)
    A = pH_ewf_2200(:,:,k,3)-pH_ews_2200(:,:,k,3);
    B = -1e-6*circ_VTna_pi(:,:,k);
    scatter(reshape(A, numel(A), 1), repmat(lambda(k), 1, numel(A)), 30, reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1), 'filled'); hold on
    ylim([.5 2.1])
    xlim([-.011 -0.0015])
    cmocean('matter')
    if k==1
       ylabel('Climate feedback parameter \lambda, W m^{-2} K^{-1}') 
       title('RCP 2.6')
    end
    
    subplot(1,3,2)
    A = pH_ewf_2200(:,:,k,1)-pH_ews_2200(:,:,k,1);
    scatter(reshape(A, numel(A), 1), repmat(lambda(k), 1, numel(A)), 30, reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1), 'filled'); hold on
    ylim([.5 2.1])
    xlim([-.011 -0.0015])
    cmocean('matter')
    if k==1
       title('RCP 4.5')
       xlabel('pH_{EW-FAST} - pH_{EW-SLOW}')
    end
    
    subplot(1,3,3)
    A = pH_ewf_2200(:,:,k,2)-pH_ews_2200(:,:,k,2);
    scatter(reshape(A, numel(A), 1), repmat(lambda(k), 1, numel(A)), 30, reshape(-1e-6*circ_VTna_pi(:,:,k), numel(A), 1), 'filled'); hold on
    ylim([.5 2.1])
    xlim([-.011 -0.0015])
    cmocean('matter')
    if k==1
        title('RCP 8.5')
    end
end

%% pH ratios
figure;
c = cmocean('-thermal', length(lambda));
for k=1:length(lambda)
    subplot(1,3,1)
    A = (pH_ews_2200(:,:,k,3) - pH_ctrl_2200(:,:,k,3)) ./ ...
        (pH_ewf_2200(:,:,k,3) - pH_ctrl_2200(:,:,k,3)) - 1;
    B = -1e-6*circ_VTna_pi(:,:,k);
    scatter(reshape(A, numel(A), 1), reshape(B, numel(B), 1), 30, c(k,:)); hold on
    if k==1
       ylabel('Northern Sinking Transport, Sv') 
       title('RCP 2.6')
    end
    
    
    subplot(1,3,2)
    A = (pH_ews_2200(:,:,k,1) - pH_ctrl_2200(:,:,k,1)) ./ ...
        (pH_ewf_2200(:,:,k,1) - pH_ctrl_2200(:,:,k,1)) - 1;
    scatter(reshape(A, numel(A), 1), reshape(B, numel(B), 1), 30, c(k,:)); hold on
    if k==1
       title({'Proportion overshoot of slow EW', 'RCP 4.5'})
       xlabel('pH_{slow EW} / pH_{fast EW}')
    end
    
    subplot(1,3,3)
    A = (pH_ews_2200(:,:,k,2) - pH_ctrl_2200(:,:,k,2))./ ...
        (pH_ewf_2200(:,:,k,2) - pH_ctrl_2200(:,:,k,2)) - 1;
    scatter(reshape(A, numel(A), 1), reshape(B, numel(B), 1), 30, c(k,:)); hold on
    if k==1
        title('RCP 8.5')
    end
end

cb = colorbar
cmocean('-thermal', length(lambda));
set(cb, 'YTick', [1/(length(lambda)*2):1/length(lambda):1], 'YTickLabel', ...
    {'0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0'});
title(cb,'\lambda (W m^{-2})')



