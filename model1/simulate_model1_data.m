% ***
% A script for simulating the data for the model 1 (performance
% of a single BCI in a group of subjects)
% ***

%% Define experiment parameters
SAVE_RESULTS = true; 
PLOT_RESULTS = true; 
FILENAME_EXP = 'SimulatedExp1'; % Name for saving simulated data (without extension)
N_S = 10; % Number of subjects in the experiment
T = repmat(100, N_S, 1); % Number of trials for each subject

%% Define true parameter values
mu_psi_true = 0.75; % Average accuracy in the population (on probability scale)
mu_alpha_true = log(mu_psi_true / (1 - mu_psi_true)); % Average accuracy in the population (on logit scale)
sigma_alpha_true = 0.5; % Accuracy SD in the population (on logit scale)

%% Generate data
% Generate subject-wise accuracies
alpha_true = normrnd(mu_alpha_true, sigma_alpha_true, N_S, 1); % Accuracy on the logit scale
psi_true = 1 ./ (1 + exp(-alpha_true)); % Accuracy on the probability scale
y = binornd(T, psi_true); % Subject-wise number of correct trials

%% Plot data
if PLOT_RESULTS
    % Prepare data to plot
    sampleAcc = y ./ T; % Compute sample accuracy
    
    % Plot figure
    figure
    hold on
    
    % Plot true subject-wise accuracies and sample accuracies
    bar([sampleAcc, psi_true])
    
    % Customize figure
    ylim([0, 1])
    legend('Sample acc.', 'True acc.')
end

%% Save data
if SAVE_RESULTS
    % Save the true parameters and data in Matlab format
    save([FILENAME_EXP '.mat'],...
        'N_S', 'T', ... % experimental design
        'mu_psi_true', 'mu_alpha_true', 'sigma_alpha_true',... % true top-level parameter values
        'psi_true', 'alpha_true',... % true subject-level parameter values
        'y'); % generated data
        
    % Save data in the format expected by 'run_model1.m' script
    dlmwrite([FILENAME_EXP '.csv'], [y T], ';', 1, 1);
end