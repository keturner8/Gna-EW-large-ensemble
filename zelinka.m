%% Script to create figure of CMIP5 and CMIP6 climate feedback parameters 
%% according to the data found in Zelinka et al. (2019). 
%% This figure is provided in the manuscript appendix

% climate feedback parameters for cmip5 vs cmip6 models
lambda_cmip5 = -1*[0.76 0.81 0.92 1.18 1.14 0.64 1.03 0.84 0.91 0.76 1.23 ...
    1.37 1.66 1.76 0.63 0.75 0.80 1.02 0.92 1.53 1.13 1.19 1.23 1.22 1.11 ...
    1.15 1.19 1.43];

lambda_cmip6 = -1*[1.03 0.92 1.82 0.63 0.71 0.74 0.92 0.62 0.65 0.63 0.81 ...
    0.78 1.40 0.82 1.45 1.13 0.63 1.48 0.75 1.54 1.40 1.22 1.10 0.76 1.34 ...
    1.05 0.67];

edges = (0.5:0.2:2.1);
%%
figure;
subplot(1,2,1)
histogram(-1* lambda_cmip5, edges)
text(1.8,10, 'N=28')
ylim([0 10.5])
ylabel('Number of models')
xlabel('\lambda (W m^{-2} K^{-1})')
t1 = title('(a) CMIP5 climate feedback parameters')
t1.Units = 'Normalize'; 
t1.Position(1) = 0; % use negative values (ie, -0.1) to move further left
t1.HorizontalAlignment = 'left';  

subplot(1,2,2)
histogram(-1* lambda_cmip6, edges)
text(1.8,10, 'N=27')
ylim([0,10.5])
ylabel('Number of models')
xlabel('\lambda (W m^{-2} K^{-1})')
t2 = title('(b) CMIP6 climate feedback parameters')
t2.Units = 'Normalize'; 
t2.Position(1) = 0; % use negative values (ie, -0.1) to move further left
t2.HorizontalAlignment = 'left';  
