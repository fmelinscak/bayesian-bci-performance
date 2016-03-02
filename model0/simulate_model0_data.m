% ***
% A script for simulating the data for the model 0 (simple
% example of normal data)
% ***

%% Define experiment parameters
SAVE_RESULTS = false; 
PLOT_RESULTS = true; 
FILENAME_EXP = 'SimulatedExp0'; % Name for saving simulated data (without extension)
N = 14; % Number of samples collected in the experiment

%% Define true parameter values
mu_true = 1; % true population mean
sigma_true = 3; % true population standard deviation

%% Generate data to illustrate difference between p-values and CIs
% The p-value should be >0.05 and CI relatively wide for illustration
while true
    y = normrnd(mu_true, sigma_true, N, 1); % generate data
    [nullRejected, p_val, ci, stats] = ttest(y);
    if ~nullRejected && (ci(2) > 3)
        break;
    end
end

% Display results for the selected dataset
nullRejected, p_val, ci, stats

%% Plot data
if PLOT_RESULTS
    % Plot figure
    figure
    hold on
    
    % Plot the data
    hist(y)
    
    hold off
end

%% Save data
if SAVE_RESULTS
    % Save the true parameters and data in Matlab format
    save([FILENAME_EXP '.mat'],...
        'N', ... % experimental design
        'mu_true', 'sigma_true',... % true parameter values
        'y'); % generated data
        
    % Save data in the format expected by 'run_model0.m' script
    dlmwrite([FILENAME_EXP '.csv'], y, ';');
end