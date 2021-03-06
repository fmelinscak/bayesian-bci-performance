% ***
% A script for performing the MCMC inference for the model 3 (comparison of
% different BCI approaches in a within-subject design)
% ***

%% Imports
addpath(fullfile('..', 'external', 'matbugs')); % add matbugs


%% Parameters of the analysis
PATH_BUGS = fullfile('C:','WinBUGS14'); % !!! IMPORTANT: set the path to the WinBUGS installation
PATH_DATA = fullfile('./data');
FILENAME_DATA = 'Brunner2011.csv'; % filename of the data to be analyzed
FILENAME_MODEL = 'winbugs-model3.txt'; % WinBUGS model used for analysis
FILENAME_MCMC = 'Brunner2011_MCMCsample.mat'; % filename for saving the MCMC sample to disk
SAVE_SAMPLE = true; % whether to save the MCMC sample to disk

viewWinBUGS = 1; % 1 - leave WinBUGS open after sampling; 0 - close WinBUGS after sampling

paramsToMonitor = ... % Parameters whose samples are recorded during MCMC
    {'beta0', 'beta1', 'eta', 'sigma_eta', 'sigma_alpha',...  % Group-level parameters
    'alpha', 'psi',...              % Individual parameters
    'eta_pred', 'mu_pred', 'alpha_pred', 'psi_pred'};      % Predicted parameters


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
