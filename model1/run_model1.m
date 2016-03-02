% ***
% A script for performing the MCMC inference for the model 1 (performance
% of a single BCI in a group of subjects)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matbugs')); % add matbugs


%% Parameters of the analysis
PATH_BUGS = fullfile('C:','WinBUGS14'); % !!! IMPORTANT: set the path to the WinBUGS installation
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'Power2010.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'winbugs-model1.txt'; % WinBUGS model used for analysis
FILENAME_MCMC = 'Power2010_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample

viewWinBUGS = 1; % 1 - leave WinBUGS open after sampling; 0 - close WinBUGS after sampling

paramsToMonitor = ...
    {'mu_alpha', 'sigma_alpha',...  % Group-level parameters
    'alpha', 'psi',...              % Individual parameters
    'alpha_pred', 'psi_pred'};      % Predicted parameters

% MCMC parameters
nChains  = 3; % Number of parallel chains
nBurnin  = 50000; % Number of burn-in samples
nSamples = 50000;  % Number of recorded samples
nThin = 1; % Thinning factor (take every nThin-th sample)

% Hyper-parameters (parameters of the prior)
M_mu_alpha = 0; % mean of the normal prior on mu_alpha
S_mu_alpha = sqrt(2); % std. dev. of the normal prior on mu_alpha

L_sigma_alpha = 0.001; % lower bound of the uniform prior on the sigma_alpha (slightly above zero for numerical reasons)
U_sigma_alpha = 10; % upper bound of the uniform prior on the sigma_alpha


%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 1); % assumes first row is the header and first column are subject labels
y = data(:, 1); % assumes the second column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 2); % assumes the third column of the CSV are the total number of trials (per subject)
N_S = size(y, 1); % Number of subjects

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;

%% Create the structure containing observed data and hyperparameters
dataStruct = struct(...
    'y', y,...
    'T', T,...
    'N_S', N_S,...
    'M_mu_alpha', M_mu_alpha,...
    'S_mu_alpha', S_mu_alpha,...
    'L_sigma_alpha', L_sigma_alpha,...
    'U_sigma_alpha', U_sigma_alpha);

%% Set initial values for the stochastic, unobserved parameters
for i = 1 : nChains
    S.alpha = repmat(0, N_S, 1);
    S.mu_alpha = normrnd(M_mu_alpha, S_mu_alpha);
    S.sigma_alpha = unifrnd(L_sigma_alpha, U_sigma_alpha);
    S.alpha_pred = normrnd(M_mu_alpha, S_mu_alpha);
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
        'samples', 'stats', 'nChains', 'nSamples')
end

