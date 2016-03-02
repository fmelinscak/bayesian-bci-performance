% ***
% A script for visualizing the results of MCMC inference for the model 1 (performance
% of a single BCI in a group of subjects)
% ***

%% Imports
addpath(fullfile('..', 'external', 'hline_vline')); % plotting vertical and horizontal lines
addpath(fullfile('..', 'external', 'kde2d')); % 2D kernel density estimation

%% Visualization parameters
PATH_DATA = fullfile('./data'); % path to analyzed data (CSV)
PATH_SAMPLE = fullfile('.'); % path to the save MCMC sample (mat)
FILENAME_DATA = 'Power2010.csv';
FILENAME_SAMPLE = 'Power2010_MCMCsample.mat';
SIMULATED_DATA = false; % If the data is simulated, inferred parameters are compared with true parameters
if SIMULATED_DATA
    PATH_TRUE_PARAMS = fullfile('./data'); % path to the true parameters (mat)
    FILENAME_TRUE_PARAMS = 'SimulatedExp1.mat';
end


%% Load data
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 1); % assumes first row is the header and first column are subject labels
y = data(:, 1); % assumes the second column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 2); % assumes the third column of the CSV are the total number of trials (per subject)
N_S = size(y, 1); % Number of subjects

if SIMULATED_DATA % Load true parameters
    load(fullfile(PATH_TRUE_PARAMS, FILENAME_TRUE_PARAMS), ...
        'mu_psi_true', 'mu_alpha_true', 'sigma_alpha_true',... % true top-level parameter values
        'psi_true', 'alpha_true'); % true subject-level parameter values
end

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;

%% Load MCMC sample
load(fullfile(PATH_SAMPLE, FILENAME_SAMPLE),...
        'samples', 'stats', 'nChains', 'nSamples');
    
%% Visulize the subjects' sample accuracies and inferred posterior accuracies
% Prepare data to plot
psi_pooled = reshape(samples.psi, nChains * nSamples, N_S); % pooling the samples accross chains
psi_prctiles = prctile(psi_pooled, [2.5 50 97.5], 1)'; % percentiles of the posterior (N_S x 3)

[psi_prctiles_srtd, iSrt] = sortrows(psi_prctiles, 2); % subjects sorted by the posterior median

iSubj = 1 : N_S;
vals = psi_prctiles_srtd(:,2); % median
err_hi = psi_prctiles_srtd(:,3) - psi_prctiles_srtd(:,2); % 97.5th - median
err_lo = psi_prctiles_srtd(:,2) - psi_prctiles_srtd(:,1); % median - 2.5th

sampleAcc_srtd = sampleAcc(iSrt);

% Plot the figure
figure
hold on

% Plot inferred subject-wise parameters
hPost = errorbar(iSubj, vals, err_lo, err_hi, 'o', 'Color', [ 0.6510    0.8078    0.8902], 'MarkerSize', 5, 'MarkerFaceColor',  [ 0.6510    0.8078    0.8902]);

% Plot sample accuracies
hData = plot(iSubj, sampleAcc_srtd, '>', 'MarkerEdgeColor','k', 'MarkerFaceColor','r','MarkerSize', 5);

% Plot true parameter values (if available)
if SIMULATED_DATA
    hTrue = plot(iSubj, psi_true(iSrt), '<', 'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize', 5);
end

% Customize the figure
ylim([0,1])
xlim([0.5, N_S+0.5])
title('Observed and inferred subject-wise acc.')
xlabel('Subjects (sorted by accuracy)')
ylabel('Accuracy')
hHline = hline(0.5, 'k--');

legendEntries = [hData hPost];
legendLabels = {'Observed sample accuracy', 'Posterior median with 95% CI'};
if SIMULATED_DATA
    legendEntries = [legendEntries hTrue];
    legendLabels = {legendLabels{:}, 'True accuracy'};
end
hLgnd = legend(legendEntries, legendLabels);
set(hLgnd, 'location', 'best', 'box', 'off')

set(gca, 'TickDir', 'out')

hold off

%% Visualize the joint posterior of mu_alpha and sigma_alpha
% Prepare data to plot
mu_alpha_pooled = reshape(samples.mu_alpha, nChains * nSamples, 1); % pooling the samples accross chains
sigma_pooled = reshape(samples.sigma_alpha, nChains * nSamples, 1); % pooling the samples accross chains

[bandwidth, pdf_density, pdf_mu_alpha, pdf_sigma_alpha] = kde2d([mu_alpha_pooled, sigma_pooled]); % 2D kernel density estimation

% Plot the figure
figure
hold on

% Plot the joint posterior
[~, hJoint] = contour3(pdf_mu_alpha, pdf_sigma_alpha, pdf_density, 10);

% Plot the true parameter values (if available)
if SIMULATED_DATA
   hTrueMean = vline(mu_alpha_true, 'k--');
   hTrueSD = hline(sigma_alpha_true, 'k--');
end
   

% Customize the figure
colormap(hot)
title('Joint posterior of \mu and \sigma')
xlabel('Group mean \mu_\alpha')
ylabel('Group std. dev. \sigma_\alpha')
set(gca, 'TickDir', 'out')
hold off


%% Visualize the posterior and the posterior predictive distributions for the group mean
% Prepare data to plot
mu_psi_pooled = 1 ./ (1 + exp(-mu_alpha_pooled)); % from the logit scale to probability scale (!! can over/under-flow for large values)
mu_psi_prctiles = prctile(mu_psi_pooled, [2.5 50 97.5], 1)'; % percentiles of the posterior


psi_pred_pooled =  reshape(samples.psi_pred, nChains * nSamples, 1); % pooling the samples accross chains
psi_pred_prctiles = prctile(psi_pred_pooled, [2.5 50 97.5], 1)';

xi = linspace(0, 1, 400);
[f_mu_psi] = ksdensity(mu_psi_pooled, xi, 'support', [-0.001, 1.001], 'function', 'pdf');
[f_psi_pred] = ksdensity(psi_pred_pooled, xi, 'support', [-0.001, 1.001], 'function', 'pdf');

maxDensity = max(max(f_mu_psi), max(f_psi_pred)); % used for scaling

% Plot the figure
figure
hold on

% Plot the data
hData = scatter(sampleAcc, 0.02 * maxDensity * ones(size(sampleAcc)), 'xr');

% Plot the posterior of the accuracy
hPost = plot(xi, f_mu_psi, 'LineWidth', 1, 'Color', 'k');
hPostCi = line(mu_psi_prctiles([1,3]), 0.08 * maxDensity * [1 1], 'Color', [0 0 0], 'LineWidth', 1.5);

% Plot the posterior predictive distribution
hPred = plot(xi, f_psi_pred, 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');
hPredCi = line(psi_pred_prctiles([1,3]), 0.05 * maxDensity * [1 1], 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle', ':');

% Plot the true mean (if available)
if SIMULATED_DATA
    hTrue = vline(mu_psi_true, 'k-.');
end

% Customize figure
set(gca, 'TickDir', 'out')
title('Posterior for \mu and posterior predictive for acc.')
xlabel('Accuracy')
xlim([0, 1])
ylabel('Probability density')
vline(0.5, 'k--')
text(0.3, 4, 'Chance level')

legendEntries = [hData hPost hPostCi hPred hPredCi];
legendLabels = {'Observed data $y_i / T_i$', '$p(\mu_\psi|d)$', '95\% posterior CI', '$p(\tilde{\psi}|d)$', '95\% predictive CI'};
if SIMULATED_DATA
    legendEntries = [legendEntries hTrue];
    legendLabels = {legendLabels{:}, 'True average accuracy'};
end
hLgnd = legend(legendEntries, legendLabels);

set(hLgnd, 'Location', 'NorthWest', 'Interpreter', 'LaTeX', 'box', 'off')

hold off
