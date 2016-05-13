% ***
% A script for performing the MCMC inference for the model 0 (simple
% example of normal data)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matbugs')); % add matbugs


%% Parameters of the analysis
PATH_BUGS = fullfile('C:','WinBUGS14'); % !!! IMPORTANT: set the path to the WinBUGS installation
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'SimulatedExp0.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'winbugs-model0.txt'; % WinBUGS model used for analysis
FILENAME_MCMC = 'SimulatedExp0_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample

viewWinBUGS = 1; % 1 - leave WinBUGS open after sampling; 0 - close WinBUGS after sampling

paramsToMonitor = ... % Parameters recorded in MCMC
    {'mu', 'sigma',...  % Population parameters
    'y_pred'};      % Predicted data

% MCMC parameters
nChains  = 3; % Number of parallel chains
nBurnin  = 50000; % Number of burn-in samples
nSamples = 50000;  % Number of recorded samples
nThin = 1; % Thinning factor (take every nThin-th sample)

% Hyper-parameters (parameters of the prior)
M_mu = 0; % mean of the normal prior on mu
S_mu = 3; % std. dev. of the normal prior on mu

L_sigma = 0.001; % lower bound of the uniform prior on the sigma_alpha (slightly above zero for numerical reasons)
U_sigma = 10; % upper bound of the uniform prior on the sigma_alpha


%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';'); % assumes data y is in the first column
y = data(:, 1); 
N = size(y, 1); % Number of samples


%% Create the structure containing observed data and hyperparameters
dataStruct = struct(...
    'y', y,...
    'N', N,...
    'M_mu', M_mu,...
    'S_mu', S_mu,...
    'L_sigma', L_sigma,...
    'U_sigma', U_sigma);


%% Set initial values for the stochastic, unobserved parameters
for i = 1 : nChains
    S.mu = normrnd(M_mu, S_mu);
    S.sigma = unifrnd(L_sigma, U_sigma);
    S.y_pred = normrnd(M_mu, S_mu);
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

