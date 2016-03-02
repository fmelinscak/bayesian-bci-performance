% ***
% A script for simulating the data for the model 2 (performance
% of a single BCI in a group of subjects)
% ***

%% Define experiment parameters
SAVE_RESULTS = true; 
PLOT_RESULTS = true; 
FILENAME_EXP = 'SimulatedExp2'; % Name for saving simulated data (without extension)
N_S = 80; % Number of subjects in the experiment
T = repmat(300, N_S, 1); % Number of trials for each subject

%% Define true parameter values
mean_x_true = 8; % Mean of the covariate
sigma_x_true = 1; % Standard deviation of the covariate
beta0_prob_true = 0.75; % Accuracy for the average value of the covariate (on probability scale)
beta0_true = log(beta0_prob_true / (1 - beta0_prob_true)); % Accuracy for the average value of the covariate (on logit scale)
beta1_true = 0.4; % Increase in logit accuracy for an increase of 1 sample SD in the covariate
sigma_alpha_true = 0.3; % Accuracy SD in the population (on logit scale)

%% Generate data
% Generate subject-wise covariate values
x = normrnd(mean_x_true, sigma_x_true, N_S, 1); % Assumption of normal distribution of the covariate
[z, m_x, s_x] = zscore(x); % Calculate the z-score, sample mean, and sample SD for the covariate

% Generate subject-wise accuracies
mu_true = beta0_true + beta1_true .* z; % predicted accuracy on the logit scale
alpha_true = normrnd(mu_true, sigma_alpha_true, N_S, 1); % accuracies on the logit scale
psi_true = 1 ./ (1 + exp(-alpha_true)); % accuracies on the probability scale
y = binornd(T, psi_true);

%% Plot data
if PLOT_RESULTS
    % Prepare data to plot
    sampleAcc = y ./ T; % Compute sample accuracy
    
    x_pred = linspace(4, 12, 200); % Values of the covariate for which to compute the true accuracy
    z_pred = zscore(x_pred);
    mu_pred = beta0_true + beta1_true .* z_pred;
    psi_pred = 1 ./ (1 + exp(-mu_pred));
    
    % Plot figure
    figure
    hold on
    
    % Plot sample accuracies on probability scale
    hData = scatter(x, sampleAcc); 
    
    % Customize figure
    ylim([0, 1])
    legend([hData], 'Sample acc.')
end

%% Save data
if SAVE_RESULTS
    % Save the true parameters and data in Matlab format
    save([FILENAME_EXP '.mat'],...
        'N_S', 'T', ... % experimental design
        'beta0_prob_true', 'beta0_true', 'beta1_true', 'sigma_alpha_true',... % true top-level parameter values
        'mu_true', 'psi_true', 'alpha_true',... % true subject-level parameter values
        'y'); % generated data
        
    % Save data in the format expected by 'run_model2.m' script
    dlmwrite([FILENAME_EXP '.csv'], [x y T], ';', 1, 1);
end