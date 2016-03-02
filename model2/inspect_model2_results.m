% ***
% A script for visualizing the results of MCMC inference for the model 2 (association 
% between a subject-specific variable and BCI performance)
% ***

%% Imports
addpath(fullfile('..', 'external', 'hline_vline')); % plotting vertical and horizontal lines
addpath(fullfile('..', 'external', 'kde2d')); % 2D kernel density estimation

%% Visualization parameters
PATH_DATA = fullfile('./data'); % path to analyzed data (CSV)
PATH_SAMPLE = fullfile('.'); % path to the save MCMC sample (mat)
FILENAME_DATA = 'Blankertz2010.csv';
FILENAME_SAMPLE = 'Blankertz2010_MCMCsample.mat';
SIMULATED_DATA = false; % If the data is simulated, inferred parameters are compared with true parameters
if SIMULATED_DATA
    PATH_TRUE_PARAMS = fullfile('./data');
    FILENAME_TRUE_PARAMS = 'SimulatedExp2.mat';
end


%% Load the original data and compute sample statistics
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 1); % assumes first row is the header and first column are subject labels
x = data(:, 1); % assumes the second column of the CSV are the subject-specific values of the covariate
y = data(:, 2); % assumes the third column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 3); % assumes the fourth column of the CSV are the total number of trials (per subject)
N_S = size(y, 1); % Number of subjects

if SIMULATED_DATA % Load true parameters
    load(fullfile(PATH_TRUE_PARAMS, FILENAME_TRUE_PARAMS), ...
        'beta0_prob_true', 'beta0_true', 'beta1_true', 'sigma_alpha_true',... % true top-level parameter values
        'mu_true', 'psi_true', 'alpha_true'); % true subject-level parameter values
end

% Compute the standardized values of the covariate
[z, m_x, s_x] = zscore(x);

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;
    

%% Load the MCMC sample and values for prediction
load(fullfile(PATH_SAMPLE, FILENAME_SAMPLE),...
        'samples', 'stats', 'nChains', 'nSamples', 'z_pred', 'N_P');


%% Visualize the data: sample accuracy vs the covariate
% Plot the figure
figure
hold on

% Plot the data
hData = scatter(z, sampleAcc, 'k');


% Customize the figure
title('Observed sample accuracies with OLS fit')
xlabel('Covariate (z-score)')
ylabel('Accuracy')
xlim([-4,4])
lsline % least square line
ylim([0 1.05])
hHline = hline(0.5, 'k--');
text(2.5, 0.52, 'Chance level')

set(gca, 'TickDir', 'out')

hold off

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

% Plot the inferred subject-wise parameters
hPost = errorbar(iSubj, vals, err_lo, err_hi, 'o', 'Color', [ 0.6510    0.8078    0.8902], 'MarkerSize', 5, 'MarkerFaceColor', [ 0.6510    0.8078    0.8902]);

% Plot sample accuracies
hData = plot(iSubj, sampleAcc_srtd, '>', 'MarkerEdgeColor','k', 'MarkerFaceColor','r','MarkerSize', 5);

% Plot true parameter values (if available)
if SIMULATED_DATA
    hTrue = plot(iSubj, psi_true(iSrt), '<', 'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize', 5);
end

% Customize the figure
title('Subject-wise inferences')
ylim([0,1])
xlim([0.5, N_S+0.5])
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


%% Visualize the joint posterior of beta0 and beta1
% Prepare data to plot
beta0_pooled = reshape(samples.beta0, nChains * nSamples, 1); % pooling the samples accross chains
beta1_pooled = reshape(samples.beta1, nChains * nSamples, 1); % pooling the samples accross chains

[bandwidth, pdf_density, pdf_beta0, pdf_beta1] = kde2d([beta0_pooled, beta1_pooled]); % 2D kernel density estimation

% Plot the figure
figure
hold on

% Plot the joint posterior
[~, hJoint] = contour3(pdf_beta0, pdf_beta1, pdf_density, 10);

% Plot the true parameter values (if available)
if SIMULATED_DATA
   hTrueIntercept = vline(beta0_true, 'k--');
   hTrueSlope = hline(beta1_true, 'k--');
end

% Customize the figure
colormap(hot)
title('Joint posterior of \beta_0 and \beta_1')
xlabel('Intercept \beta_0')
ylabel('Slope \beta_1')
set(gca, 'TickDir', 'out')
hold off

%% Visualize the marginal posteriors of the intercept and the slope
% Prepare data to plot
beta0_pooled = reshape(samples.beta0, nChains * nSamples, 1); % pooling the samples accross chains
beta0_prctiles = prctile(beta0_pooled, [2.5 50 97.5], 1)';

beta1_pooled =  reshape(samples.beta1, nChains * nSamples, 1); % pooling the samples accross chains
beta1_prctiles = prctile(beta1_pooled, [2.5 50 97.5], 1)';

xi = linspace(-5, 5, 400); % points for which the density is estimated
[f_beta0] = ksdensity(beta0_pooled, xi, 'support', 'unbounded', 'function', 'pdf', 'width', 0.04);
[f_beta1] = ksdensity(beta1_pooled, xi, 'support', 'unbounded', 'function', 'pdf', 'width', 0.04);


max_f_beta0 = max(f_beta0);
max_f_beta1 = max(f_beta1);

% Plot the figure
figure
hold on

% Plot the posterior of the intercept
subplot(2,1,1)
hPost_beta0 = plot(xi, f_beta0, 'LineWidth', 1, 'Color', 'k');
hPost_beta0_ci = line(beta0_prctiles([1,3]), max_f_beta0 * [0.05 0.05], 'Color', [0 0 0], 'LineWidth', 1.5);

% Plot the true intercept value (if available)
if SIMULATED_DATA
    hTrueIntercept = vline(beta0_true, 'r:');
end

% Customize the subplot
title('Marginal posterior for \beta_0')
xlabel('Intercept \beta_0 (logit scale)')
ylabel('p(\beta_0 | d)')
hVline = vline(0, 'k--');
text(-0.5, max_f_beta0 * 0.2, 'Chance level', 'Rotation', 90)
xlim([-5, 5])
ylim([0, max_f_beta0 * 1.3])
set(gca, 'TickDir', 'out')

% Plot the posterior of the slope
subplot(2,1,2)
hPost_beta1 = plot(xi, f_beta1, 'LineWidth', 1, 'Color', 'k');
hPost_beta1_ci = line(beta1_prctiles([1,3]), max_f_beta1 * [0.05 0.05], 'Color', [0 0 0], 'LineWidth', 1.5);

% Plot the true intercept value (if available)
if SIMULATED_DATA
    hTrueSlope = vline(beta1_true, 'r:');
end

% Customize the subplot
title('Marginal posterior for \beta_1')
xlabel('Slope \beta_1 (logit scale)')
ylabel('p(\beta_1 | d)')
hVline = vline(0, 'k--');
text(-2.5, max_f_beta1 * 1.2, 'Negative effect <')
text(0.25, max_f_beta1 * 1.2, '> Positive effect')
xlim([-5, 5])
ylim([0, max_f_beta1 * 1.3])

% Customize figure

set(gca, 'TickDir', 'out')

hold off

%% Visualize the data and the posterior predictive distribution together
% Prepare data to plot
mu_pred_pooled = reshape(samples.mu_pred, nChains * nSamples, N_P); % pooling the samples accross chains
mu_pred_prob_pooled = 1 ./ (1 + exp(-mu_pred_pooled)); % from the logit scale to probability scale (!! can over/under-flow for large values)
mu_pred_prob_prctiles = prctile(mu_pred_prob_pooled, [2.5 50 97.5], 1)';

psi_pred_pooled = reshape(samples.psi_pred, nChains * nSamples, N_P);
psi_pred_prctiles = prctile(psi_pred_pooled, [2.5 50 97.5], 1)';

if SIMULATED_DATA
    psi_pred_true = 1 ./ (1 + exp(-(beta0_true + beta1_true * z_pred)));
end

% Plot the figure
figure
hold on

% Plot the data
hData = scatter(z, sampleAcc, 'k');

% Plot the posterior median and 95% CI
hPosteriorMedian = plot(z_pred, mu_pred_prob_prctiles(:,2),...
    'Color', [0.9882    0.5529    0.3843], 'LineWidth', 1.5);
hPosterior95ci = plot(z_pred, mu_pred_prob_prctiles(:,[1 3]),...
    'Color', [0.9882    0.5529    0.3843], 'LineWidth', 1.5, 'LineStyle', '--');

% Plot the posterior predictive 95% CI
hPredictive95ci = plot(z_pred, psi_pred_prctiles(:,[1 3]),...
    'Color', [0.5529    0.6275    0.7961], 'LineWidth', 1.5, 'LineStyle', ':');

% Plot the true mean (if available)
if SIMULATED_DATA
    hTrueMean = plot(z_pred, psi_pred_true,...
        'Color', [0.4000    0.7608    0.6471], 'LineWidth', 1.5);
end

% Customize figure
xlabel('Covariate (z-score)')
ylabel('Accuracy')
ylim([0 1.05])
xlim([-4,4])
hHline = hline(0.5, 'k--');
text(2.5, 0.52, 'Chance level')

set(gca, 'TickDir', 'out')

hold off

