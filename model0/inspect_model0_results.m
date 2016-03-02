% ***
% A script for visualizing the results of MCMC inference for the model 0 (simple
% example of normal data)
% ***

%% Imports
addpath(fullfile('..', 'external', 'hline_vline')); % plotting vertical and horizontal lines
addpath(fullfile('..', 'external', 'kde2d'));  % 2D kernel density estimation

%% Visualization parameters
PATH_DATA = fullfile('./data'); % path to analyzed data (CSV)
PATH_SAMPLE = fullfile('.'); % path to the save MCMC sample (mat)
FILENAME_DATA = 'SimulatedExp0.csv';
FILENAME_SAMPLE = 'SimulatedExp0_MCMCsample.mat';
SIMULATED_DATA = true; % If the data is simulated, inferred parameters are compared with true parameters
if SIMULATED_DATA
    PATH_TRUE_PARAMS = fullfile('./data'); % path to the true parameters (mat)
    FILENAME_TRUE_PARAMS = 'SimulatedExp0.mat';
end


%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';'); % assumes data y is in the first column
y = data(:, 1); 
N = size(y, 1); % Number of samples

if SIMULATED_DATA % Load true parameters
    load(fullfile(PATH_TRUE_PARAMS, FILENAME_TRUE_PARAMS), ...
        'mu_true', 'sigma_true'); % true population parameter values     
end


%% Load MCMC sample
load(fullfile(PATH_SAMPLE, FILENAME_SAMPLE),...
        'samples', 'stats', 'nChains', 'nSamples');
    

%% Visualize the joint posterior of mu and sigma
% Prepare data to plot
mu_pooled = reshape(samples.mu, nChains * nSamples, 1); % pooling the samples accross chains
sigma_pooled = reshape(samples.sigma, nChains * nSamples, 1); % pooling the samples accross chains

[bandwidth, pdf_density, pdf_mu, pdf_sigma] = kde2d([mu_pooled, sigma_pooled]); % 2D kernel density estimation

% Plot the figure
figure
hold on

% Plot the joint posterior
[~, hJoint] = contour3(pdf_mu, pdf_sigma, pdf_density, 10);

   
% Customize the figure
colormap(hot)
title('Full posterior distribution')
xlabel('Population mean \mu')
ylabel('Population std. dev. \sigma')
set(gca, 'TickDir', 'out')

% Plot the true parameter values (if available)
if SIMULATED_DATA
   hTrueMean = vline(mu_true, 'k--');
   hTrueSD = hline(sigma_true, 'k--');
end

hold off

%% Visualize the marginal posterior and prior of population mean mu
% Define the hyper-parameters of the normal prior that were used for the analysis
M_mu = 0;
S_mu = 3;

% Prepare data to plot
mu_pooled = reshape(samples.mu, nChains * nSamples, 1); % pooling the samples accross chains

mu_prctiles = prctile(mu_pooled, [2.5 50 97.5], 1)'; % Compute the median and the 95% CI

xi = linspace(min(min(mu_pooled), -3*S_mu + M_mu), max(max(mu_pooled), 3*S_mu + M_mu), 400);
[f_mu] = ksdensity(mu_pooled, xi, 'function', 'pdf', 'width', 0.4); % Marginal density estimation



% Evaluate prior at points xi
f_mu_prior = normpdf(xi, M_mu, S_mu);

% Plot the figure
figure
hold on

% Plot the marginal posterior of mu
hPost = plot(xi, f_mu, 'LineWidth', 1, 'Color', 'k');
hPostCi = line(mu_prctiles([1,3]), max(f_mu) * 0.05 * [1 1], 'Color', [0 0 0], 'LineWidth', 1.5);

% Plot the marginal prior of mu
hPrior = plot(xi, f_mu_prior, 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');

% Plot the true mean (if available)
if SIMULATED_DATA
    hTrue = vline(mu_true, 'k-.');
end

% Customize figure
set(gca, 'TickDir', 'out')
title('Marginal prior and posterior of mean \mu')
xlabel('Population mean \mu')
xlim([min(xi), max(xi)])
ylabel('Probability density')


legendEntries = [hPrior hPost hPostCi];
legendLabels = {'Prior $p(\mu)$', 'Posterior $p(\mu|d)$', '95\% posterior CI'};
if SIMULATED_DATA
    legendEntries = [legendEntries hTrue];
    legendLabels = {legendLabels{:}, 'True population mean'};
end
hLgnd = legend(legendEntries, legendLabels);

set(hLgnd, 'Location', 'NorthWest',  'Interpreter', 'LaTeX', 'box', 'off')

hold off


%% Visualize the posterior and the posterior predictive distributions for the group mean
% Prepare data to plot
[yFreq, yi] = hist(y); % The data histogram

y_pred_pooled =  reshape(samples.y_pred, nChains * nSamples, 1); % pooling the samples accross chains
y_pred_prctiles = prctile(y_pred_pooled, [2.5 50 97.5], 1)';

xi = linspace(min(min(y_pred_pooled), min(yi)), max(max(y_pred_pooled), max(yi)));
[f_y_pred] = ksdensity(y_pred_pooled, xi, 'function', 'pdf', 'width', 0.4); % Posterior predictive density estimation

densityScaling = max(yFreq) / max(f_y_pred); % scale the posterior predictive distribution to be commensurate with the histogram (visual adjustment)

% Plot the figure
figure
hold on

% Plot the histogram of the data
hHist = bar(yi, yFreq, 'hist');


% Plot the data
hData = scatter(y, max(yFreq)* 0.05 * ones(size(y)), 'xr');


% Plot the posterior predictive distribution
hPred = plot(xi, densityScaling .* f_y_pred, 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');
% hPredCi = line(y_pred_prctiles([1,3]), [0.7 0.7], 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle', ':');


% Customize figure
set(hHist, 'FaceColor', [0.85 0.85 0.85])
set(gca, 'TickDir', 'out', 'YTick', [0, 1, 2, 3])
title('Observed data and the posterior predictive distribution')
xlabel('y')
xlim([min(xi), max(xi)])
ylim([0, max(yFreq) * 1.3])
ylabel('Frequency')

legendEntries = [hData hPred];
legendLabels = {'Observed data $y_i$', 'Posterior predictive distribution $p(\tilde{y}|d)$'};

hLgnd = legend(legendEntries, legendLabels);

set(hLgnd, 'Location', 'NorthWest', 'Interpreter', 'LaTeX', 'box', 'off')


hold off
