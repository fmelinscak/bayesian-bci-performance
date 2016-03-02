% ***
% A script for simulating the data for the model 3 (comparison of
% different BCI approaches in a within-subject design)
% ***

%% Define experiment parameters
SAVE_RESULTS = true; 
PLOT_RESULTS = true; 
FILENAME_EXP = 'SimulatedExp3'; % Name for saving simulated data (without extension)
N_S = 14; % Number of subjects in the experiment
N_L = 3; % Number of factor levels in the experiment
N_O = N_S * N_L; % Total number of observations in the experiment
T = repmat(100, N_O, 1); % Number of trials for each subject

%% Define true parameter values
beta0_prob_true = 0.75; % Grand average accuracy (over factor levels and subjects; on probability scale)
beta0_true = log(beta0_prob_true / (1 - beta0_prob_true)); % Accuracy for the average value of the covariate (on logit scale)
beta1_true = [-0.3 0.1 0.2]; % Deviations from the grand average for each factor level
beta1_true = beta1_true - mean(beta1_true); % Ensure that the deviations satisfy the sum-to-zero constraint
sigma_alpha_true = 0.1; % SD of unexplained variation (on logit scale)
sigma_eta_true = 0.7; % SD of inter-subject variation

%% Generate data
% Make the observation labels (subject and factor level indicators)
s = reshape(repmat(1 : N_S, N_L, 1), [], 1); % subject indicator for each observation
l = repmat([1 : N_L]', N_S, 1); % level indicator for each observation

% Generate subject-wise random effects
eta_unadjusted = normrnd(0, sigma_eta_true, N_S, 1);
eta_true = eta_unadjusted - mean(eta_unadjusted); % Constrain to sum to zero

% Generate subject-wise accuracies
mu_true = NaN(N_O, 1);
alpha_true = NaN(N_O, 1);
psi_true = NaN(N_O, 1);
y = NaN(N_O, 1);

for iObs = 1 : N_O
    mu_true(iObs) = beta0_true + beta1_true(l(iObs)) + eta_true(s(iObs)); % predicted accuracy on the logit scale
    alpha_true(iObs) = normrnd(mu_true(iObs), sigma_alpha_true); % accuracies on the logit scale
    psi_true(iObs) = 1 ./ (1 + exp(-alpha_true(iObs))); % accuracies on the probability scale
    y(iObs) = binornd(T(iObs), psi_true(iObs));
end


%% Plot data
if PLOT_RESULTS
    sampleAcc = y ./ T;
    barinput = reshape(sampleAcc, N_L, N_S);
    figure, bar(barinput')
    xlabel('Subject')
    ylabel('Accuracy')
    legend(arrayfun(@(x) sprintf('Level %d', x), 1 : N_L, 'UniformOutput', false), 'Location', 'best')
end
%% Save data
if SAVE_RESULTS
    % Save the true parameters and data in Matlab format
    save([FILENAME_EXP '.mat'],...
        'N_S', 'N_L', 'N_O', 'T', 's', 'l', ... % experimental design
        'beta0_prob_true', 'beta0_true', 'beta1_true', 'sigma_alpha_true', 'sigma_eta_true',... % true top-level parameter values
        'eta_true', 'mu_true', 'psi_true', 'alpha_true',... % true subject-level parameter values
        'y'); % generated data
        
    % Save data in the format expected by 'run_model3.m' script
    dlmwrite([FILENAME_EXP '.csv'], [s l y T], ';', 1, 2);
end