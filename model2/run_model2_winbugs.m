% ***
% A script for performing the MCMC inference for the model 2 (association 
% between a subject-specific variable and BCI performance)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matbugs')); % add matbugs


%% Parameters of the analysis
PATH_BUGS = fullfile('C:','WinBUGS14'); % !!! IMPORTANT: set the path to the WinBUGS installation
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'Blankertz2010.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'winbugs-model2.txt'; % WinBUGS model used for analysis
FILENAME_MCMC = 'Blankertz2010_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample to disk

viewWinBUGS = 1; % 1 - leave WinBUGS open after sampling; 0 - close WinBUGS after sampling

paramsToMonitor = ... % Parameters whose samples are recorded during MCMC
    {'beta0', 'beta1', 'sigma_alpha',...  % Group-level parameters
    'alpha', 'psi',...              % Individual parameters
    'mu_pred', 'alpha_pred', 'psi_pred'};      % Predicted parameters


% MCMC parameters
nChains  = 3; % Number of parallel chains
nBurnin  = 50000; % Number of burn-in samples
nSamples = 50000;  % Number of recorded samples
nThin = 1; % Thinning factor (take every nThin-th sample)


% Hyper-parameters (parameters of the prior)
M_beta0 = 0; % mean of the normal prior on beta0
S_beta0 = sqrt(2); % std. dev. of the normal prior on beta0

M_beta1 = 0; % mean of the normal prior on beta1
S_beta1 = 5; % std. dev. of the normal prior on beta1

L_sigma_alpha = 0.001; % lower bound of the uniform prior on the sigma_alpha (slightly above zero for numerical reasons)
U_sigma_alpha = 10; % upper bound of the uniform prior on the sigma_alpha

% Parameters of the prediction
PREDICT_ORIGINAL_SCALE = false; % Predictions on the original scale of the covariate or on the normalized (z-score) scale
if PREDICT_ORIGINAL_SCALE
    x_pred = linspace(0, 20, 40)'; % points for prediction
    N_P = size(x_pred, 1); % number of predictions
else
    z_pred = linspace(-4, 4, 40)'; % points for prediction
    N_P = size(z_pred, 1); % number of predictions
end

%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 1); % assumes first row is the header and first column are subject labels
x = data(:, 1); % assumes the second column of the CSV are the subject-specific values of the covariate
y = data(:, 2); % assumes the third column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 3); % assumes the fourth column of the CSV are the total number of trials (per subject)
N_S = size(y, 1); % Number of subjects

% Compute the standardized values of the covariate
[z, m_x, s_x] = zscore(x);

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;

% Compute for which values predictions are made
if PREDICT_ORIGINAL_SCALE
    z_pred = (x_pred - m_x) ./ s_x;
else
    x_pred = (z_pred .* s_x) + m_x; 
end

%% Create the structure containing observed data and hyperparameters
dataStruct = struct(...
    'y', y,...
    'T', T,...
    'z', z,...
    'N_S', N_S,...
    'M_beta0', M_beta0,...
    'S_beta0', S_beta0,...
    'M_beta1', M_beta1,...
    'S_beta1', S_beta1,...
    'L_sigma_alpha', L_sigma_alpha,...
    'U_sigma_alpha', U_sigma_alpha,...
    'z_pred', z_pred, ...
    'N_P', N_P);

%% Set initial values for the stochastic, unobserved parameters
for i = 1 : nChains
    S.alpha = repmat(0, N_S, 1);
    S.beta0 = normrnd(M_beta0, S_beta0);
    S.beta1 = normrnd(M_beta1, S_beta1);
    S.sigma_alpha = unifrnd(L_sigma_alpha, U_sigma_alpha);
    S.alpha_pred = repmat(0, N_P, 1);
    initStructs(i) = S; % initStructs has the initial values for all the Markov chains
end
    
%% Calling BUGS to obtain the MCMC sample
fprintf( 'Running BUGS...\n' );
tic
[samples, stats] = matbugs(dataStruct, ...
  fullfile(pwd, FILENAME_MODEL), ...
  'init', initStructs, ...
  'nchains', nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', nSamples,...
  'nthin', nThin,...
  'view', viewWinBUGS,...
  'monitorParams', paramsToMonitor, ...
  'Bugdir', PATH_BUGS);
toc


%% Save the MCMC sample
if SAVE_SAMPLE
    save(FILENAME_MCMC,...
        'samples', 'stats', 'nChains', 'nSamples', 'z_pred', 'x_pred', 'N_P')
end

