% ***
% A script for performing the MCMC inference for the model 0 (simple
% example of normal data)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matjags')); % add matjags


%% Parameters of the analysis
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'SimulatedExp0.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'jags-model0.txt'; % JAGS model used for analysis
FILENAME_MCMC = 'SimulatedExp0_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample


paramsToMonitor = ... % Parameters recorded in MCMC
    {'mu', 'sigma',...  % Population parameters
    'y.pred'};      % Predicted data

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
    

%% Calling JAGS to obtain the MCMC sample
fprintf( 'Running JAGS...\n' );
tic
[samples, stats, ~] = matjags( ... 
    dataStruct, ...                     % Observed data
    fullfile(pwd, FILENAME_MODEL), ...    % File that contains model definition
    initStructs, ...                          % Initial values for latent variables
    'doparallel' , 0, ...      % Parallelization flag
    'nchains', nChains,...              % Number of MCMC chains
    'nburnin', nBurnin,...              % Number of burnin steps
    'nsamples', nSamples, ...           % Number of samples to extract
    'thin', nThin, ...                      % Thinning parameter
    'monitorparams', paramsToMonitor, ...     % List of latent variables to monitor
    'savejagsoutput' , 0, ...          % Save command line output produced by JAGS?
    'workingdir' , 'tmpjags' ,...
    'verbosity' , 1, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 0);                    % clean up of temporary files?
toc


%% Save the MCMC sample
if SAVE_SAMPLE
    save(FILENAME_MCMC,...
        'samples', 'stats', 'nChains', 'nSamples')
end

