% ***
% A script for performing the MCMC inference for the model 3 (comparison of
% different BCI approaches in a within-subject design)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matjags')); % add matjags


%% Parameters of the analysis
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'Brunner2011.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'jags-model3.txt'; % JAGS model used for analysis
FILENAME_MCMC = 'Brunner2011_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample to disk


paramsToMonitor = ... % Parameters whose samples are recorded during MCMC
    {'beta0', 'beta1', 'eta', 'sigma.eta', 'sigma.alpha',...  % Group-level parameters
    'alpha', 'psi',...              % Individual parameters
    'eta.pred', 'mu.pred', 'alpha.pred', 'psi.pred'};      % Predicted parameters


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

L_sigma_eta = 0.001; % lower bound of the uniform prior on the sigma_eta (slightly above zero for numerical reasons)
U_sigma_eta = 10; % upper bound of the uniform prior on the sigma_eta


%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 2); % assumes first row is the header and first two columns are subject and factor labels
s = data(:, 1); % assumes the third column of the CSV are the subject indices
l = data(:, 2); % assumes the fourth column of the CSV are the factor levels
y = data(:, 3); % assumes the fifth column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 4); % assumes the sixth column of the CSV are the total number of trials (per subject)

N_O = size(y, 1); % Total number of observations
N_S = max(s); % Number of subjects
N_L = max(l); % Number of factor levels

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;


%% Create the structure containing observed data and hyperparameters
dataStruct = struct(...
    'y', y,...
    'T', T,...
    's', s,...
    'l', l,...
    'N_O', N_O,...
    'N_S', N_S,...
    'N_L', N_L,...
    'M_beta0', M_beta0,...
    'S_beta0', S_beta0,...
    'M_beta1', M_beta1,...
    'S_beta1', S_beta1,...
    'L_sigma_alpha', L_sigma_alpha,...
    'U_sigma_alpha', U_sigma_alpha,...
    'L_sigma_eta', L_sigma_eta,...
    'U_sigma_eta', U_sigma_eta);

%% Set initial values for the stochastic, unobserved parameters
for i = 1 : nChains
    S.alpha = repmat(0, N_O, 1);
    S.beta0 = normrnd(M_beta0, S_beta0);
    S.beta1 = [NaN; normrnd(M_beta1, S_beta1, N_L - 1, 1)]; % By the STZ constraint beta1[1] is not stochastic
    S.eta = [NaN; normrnd(0, 1, N_S - 1, 1)];  % By the STZ constraint eta[1] is not stochastic
    S.sigma_eta = unifrnd(L_sigma_eta, U_sigma_eta);
    S.sigma_alpha = unifrnd(L_sigma_alpha, U_sigma_alpha);
    S.eta_pred = 0;
    S.alpha_pred = repmat(0, N_L, 1);
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
