% ***
% A script for visualizing the results of MCMC inference for the model 3 ((comparison of
% different BCI approaches in a within-subject design)
% ***

%% Imports
addpath(fullfile('..', 'external', 'hline_vline')); % plotting vertical and horizontal lines
addpath(fullfile('..', 'external', 'distributionPlot')); % distribution plots

%% Visualization parameters
PATH_DATA = fullfile('./data'); % path to analyzed data (CSV)
PATH_SAMPLE = fullfile('.'); % path to the save MCMC sample (mat)
FILENAME_DATA = 'Brunner2011.csv';
FILENAME_SAMPLE = 'Brunner2011_MCMCsample.mat';
SIMULATED_DATA = false; % If the data is simulated, inferred parameters are compared with true parameters
if SIMULATED_DATA
    PATH_TRUE_PARAMS = fullfile('./data');
    FILENAME_TRUE_PARAMS = 'SimulatedExp3.mat';
end

FACTOR_LABELS = {'ERD', 'SSVEP', 'Hybrid'}; % Labels for factor levels (for plotting purposes)


%% Load the original data and compute sample statistics
data = dlmread(fullfile(PATH_DATA, FILENAME_DATA), ';', 1, 2); % assumes first row is the header and first two columns are subject and factor labels
s = data(:, 1); % assumes the third column of the CSV are the subject indices
l = data(:, 2); % assumes the fourth column of the CSV are the factor levels
y = data(:, 3); % assumes the fifth column of the CSV are the numbers of succesful trials (per subject)
T = data(:, 4); % assumes the sixth column of the CSV are the total number of trials (per subject)

N_O = size(y, 1); % Total number of observations
N_S = max(s); % Number of subjects
N_L = max(l); % Number of factor levels

if SIMULATED_DATA % Load true parameters
    load(fullfile(PATH_TRUE_PARAMS, FILENAME_TRUE_PARAMS), ...
        'beta0_prob_true', 'beta0_true', 'beta1_true', 'sigma_alpha_true', 'sigma_eta_true',... % true top-level parameter values
        'eta_true', 'mu_true', 'psi_true', 'alpha_true'); % true subject-level parameter values
end

% Compute the sample accuracy (for later comparison with inferred accuracies)
sampleAcc = y ./ T;
    

%% Load the MCMC sample
load(fullfile(PATH_SAMPLE, FILENAME_SAMPLE),...
        'samples', 'stats', 'nChains', 'nSamples');
    
%% Visualize the data: sample accuracy vs the factor level
% Prepare the data for plotting
sampleAccTbl = NaN(N_S, N_L);

for i = 1 : N_O
    sampleAccTbl(s(i), l(i)) = sampleAcc(i); % Reorganize observations into a tabular N_S x N_L form
    
end

if SIMULATED_DATA
    psi_trueTbl = NaN(N_S, N_L);
    for i = 1 : N_O
        psi_trueTbl(s(i), l(i)) = psi_true(i); % Reorganize observations into a tabular N_S x N_L form
    end
end

% Plot the figure
figure
hold on

% Plot the data
hData = plot(1:N_L, sampleAccTbl, 'o-', 'Color', [0.5 0.5 0.5]);

% Plot the true accuracies (if available)
if SIMULATED_DATA
    hTruePsi = plot(1:N_L, psi_trueTbl, 'r+-');
end

% Customize the figure
xlabel('Factor level')
ylabel('Accuracy')
xlim([0.75 3.25])
ylim([0 1.0])
hHline = hline(0.5, 'k--', 'Chance level');

set(gca, 'XTick', 1 : N_L, 'XTickLabel', FACTOR_LABELS);

hold off

%% Visulize the subjects' sample accuracies and inferred posterior accuracies for all factor levels
% Prepare data to plot
psi_pooled = reshape(samples.psi, nChains * nSamples, N_O); % pooling the samples accross chains
psi_prctiles = prctile(psi_pooled, [2.5 50 97.5], 1)'; % percentiles of the posterior (N_O x 3)

psi_prctiles_srtd = psi_prctiles; % TODO: sorting
iSrt = 1 : N_O; % TODO: sorting
% [psi_prctiles_srtd, iSrt] = sortrows(psi_prctiles, 2); % subjects sorted by the posterior median

iObs = 1 : N_O;
vals = psi_prctiles_srtd(:,2); % median
err_hi = psi_prctiles_srtd(:,3) - psi_prctiles_srtd(:,2); % 97.5th - median
err_lo = psi_prctiles_srtd(:,2) - psi_prctiles_srtd(:,1); % median - 2.5th

sampleAcc_srtd = sampleAcc(iSrt);

% Plot the figure
figure
hold on


% Plot the posterior median and 95% CI
% hPost = errorbar(iObs, vals, err_lo, err_hi, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

hPost = errorbar(reshape(iObs, N_L, N_S)', reshape(vals, N_S, N_L), reshape(err_lo, N_S, N_L), reshape(err_hi, N_S, N_L), 'o', 'MarkerSize', 5);

% Plot the data
hData = plot(reshape(iObs, N_L, N_S)', reshape(sampleAcc_srtd, N_S, N_L), '>', 'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize', 5);

% Plot the true parameter values if available
if SIMULATED_DATA
    hTrue = plot(reshape(iObs, N_L, N_S)', reshape(psi_true, N_S, N_L), '<', 'MarkerEdgeColor','k', 'MarkerFaceColor','b','MarkerSize', 5);
end

% Customize the figure
ylim([0,1])
xlim([0.5,N_O+0.5])
xlabel('Subject')
ylabel('Accuracy')

hHline = hline(0.5, 'k--');

legendEntries = [hData(1) hPost];
condLabels = cellfun(@(s) sprintf('Posterior median with 95%% CI (cond. %s)', s), FACTOR_LABELS, 'UniformOutput', false);
legendLabels = {'Observed sample accuracy', condLabels{:}};
if SIMULATED_DATA
    legendEntries = [legendEntries hTrue(1)];
    legendLabels = {legendLabels{:}, 'True accuracy'};
end
hLgnd = legend(legendEntries, legendLabels);
set(hLgnd, 'location', 'best', 'box', 'off')

set(gca, 'TickDir', 'out', 'XTick', [1 : N_L : N_O] + floor(N_L/2), 'XTickLabel', 1 : N_S)

hold off



%% Visualize the marginal posteriors of the factor level means
% Prepare data to plot
beta0_pooled = reshape(samples.beta0, nChains * nSamples, 1); % pooling the samples accross chains
beta0_prctiles = prctile(beta0_pooled, [2.5 50 97.5], 1)';

beta1_pooled =  reshape(samples.beta1, nChains * nSamples, N_L); % pooling the samples accross chains
beta1_prctiles = prctile(beta1_pooled, [2.5 50 97.5], 1)';

level_alpha_pooled = bsxfun(@plus, beta0_pooled, beta1_pooled);

level_alpha_prctiles = prctile(level_alpha_pooled, [2.5 50 97.5], 1)';

level_psi_pooled = 1 ./ (1 + exp(-(level_alpha_pooled)));
level_psi_prctiles = prctile(level_psi_pooled, [2.5 50 97.5], 1)';

if SIMULATED_DATA
    level_alpha_true = bsxfun(@plus, beta0_true, beta1_true);
    level_psi_true = 1 ./ (1 + exp(-(level_alpha_true)));
end


% Plot the figure
figure
hold on

% Plot the level mean posteriors as violin plots
hPostViolin_level_psi = distributionPlot(level_psi_pooled,...
    'histOpt', 1,...
    'divFactor', 45,...
    'color', [0.5529    0.6275    0.7961],...
    'showMM', 0,...
    'distWidth', 0.8,...
    'xNames', FACTOR_LABELS);

% Plot the true level mean (if available)
if SIMULATED_DATA
    hTrue = plot(1 : N_L, level_psi_true, 'r.');
end

% Plot the posterior median and CI for the level mean
hPostErr_level_psi = errorbar(1 : N_L,...
    level_psi_prctiles(:,2),... % Value
    level_psi_prctiles(:,2) - level_psi_prctiles(:,1),... % Err. hi.
    level_psi_prctiles(:,3) - level_psi_prctiles(:,2), 'xk'); % Err. lo.


% Customize figure
xlim([0.5 3.5])
ylim([0 1.00])
ylabel('Accuracy: logit^{-1}(\beta_0 + \beta_{1,k})')
xlabel('Factor level')

hHline = hline(0.5, 'k--', 'Chance level');

set(gca, 'TickDir', 'out')

hold off

%% Plot the posteriors of the contrasts
% Prepare data to plot
beta1_pooled =  reshape(samples.beta1, nChains * nSamples, N_L); % pooling the samples accross chains

cntrsts = ...
    [-1 0 1; % Hybrid vs. ERD
    0 -1 1; % Hybrid vs. SSVEP
    -1 1 0; % SSVEP vs. ERD
    -0.5 -0.5 1]; % Hybrid vs. NonHybrid
cntrstLabels = {'Hybrid vs. ERD', 'Hybrid vs. SSVEP', 'SSVEP vs. ERD', 'Hybrid vs. Non-hybrid'};

nCntrsts = size(cntrsts, 1);
cntrst_vals = beta1_pooled * cntrsts';  % Compute contrsts (nChains * nSamples) x nCntrst
cntrst_prctiles = prctile(cntrst_vals, [2.5 50 97.5], 1)';

if SIMULATED_DATA
    cntrst_true = beta1_true * cntrsts';  % Compute true contrast values
end

% Plot the figure
figure
hold on

% Plot each contrast separately
xi = linspace(-10, 10, 400);
for iCntrst = 1 : nCntrsts
    subplot(2, ceil(nCntrsts / 2), iCntrst)


    f_cntrst = ksdensity(cntrst_vals(:, iCntrst), xi, 'support', 'unbounded', 'function', 'pdf', 'width', 0.1);
    
    % Plot the posterior distribution and the 
    hPost_cntrst = plot(xi, f_cntrst, 'LineWidth', 1, 'Color', 'k');
    hPost_cntrst_ci = line(cntrst_prctiles(iCntrst, [1,3]), 0.05*max(f_cntrst)*[1 1], 'Color', [0 0 0], 'LineWidth', 1.5);
    
    % Plot the true contrast value (if available)
    if SIMULATED_DATA
        hTrueContrast = vline(cntrst_true(iCntrst), 'r:');
    end
    
    % Customize subplot
    title(cntrstLabels{iCntrst})
    ylabel('p(contrast)')
    xlabel('Contrast value')

%     text(-0.5, 2.5, 'Negative effect <')
%     text(0.05, 2.5, '> Positive effect')
    xlim([-5,5])
    ylim([0, 1.1*max(f_cntrst)])
    set(gca, 'TickDir', 'out')
    h_vline = vline(0, 'k--');
    
    
%     ylim([0,5])
    
end

%% Visualize the marginal posteriors of subject's random effects
% Prepare data to plot
eta_pooled = reshape(samples.eta, nChains * nSamples, N_S); % pooling the samples accross chains
eta_prctiles = prctile(eta_pooled, [2.5 50 97.5], 1)';

vals = eta_prctiles(:,2); % median
err_hi = eta_prctiles(:,3) - eta_prctiles(:,2); % 97.5th - median
err_lo = eta_prctiles(:,2) - eta_prctiles(:,1); % median - 2.5th


% Plot the figure
figure
hold on

% Plot the posterior median and 95% CI
hPost = errorbar(1 : N_S, vals, err_lo, err_hi, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% Plot the true parameter values (if available)
if SIMULATED_DATA
    hTrue = plot(1 : N_S, eta_true, '<', 'MarkerEdgeColor','k', 'MarkerFaceColor','b','MarkerSize', 5);
end

% Customize the figure
xlim([0.5,N_S+0.5])
xlabel('Subjects (unsorted)')
set(gca, 'XTick', 1 : N_S)
ylabel('Subject random effect \eta_j (logit scale)')

hHline = hline(0, 'k--');

legendEntries = [hPost];
legendLabels = {'Posterior median with 95% CI'};
if SIMULATED_DATA
    legendEntries = [legendEntries hTrue];
    legendLabels = {legendLabels{:}, 'True value'};
end
hLgnd = legend(legendEntries, legendLabels);
set(hLgnd, 'location', 'best', 'box', 'off')

set(gca, 'TickDir', 'out')

hold off

%% Visualize the data and the posterior predictive distribution together
% Prepare data to plot
level_psi_pred_pooled = reshape(samples.psi_pred, nChains * nSamples, N_L); % pooling the samples accross chains
level_psi_pred_prctiles = prctile(level_psi_pred_pooled, [2.5 50 97.5], 1)';

yi = linspace(0, 1.001, 400);
f_level_psi_pred = zeros(length(yi), N_L);
for k = 1 : N_L
    f_level_psi_pred(:,k) = ksdensity(level_psi_pred_pooled(:,k), yi, 'support', [0, 1.001], 'function', 'pdf', 'width', 0.04);
end

% Plot the figure
figure
hold on

% Plot the data
hData = plot(1:N_L, sampleAccTbl, 'o-', 'Color', [0.5 0.5 0.5]);


% Plot the predictive distributions for all the levels
for k = 1 : N_L
    plot(k + 0.8 * f_level_psi_pred(:,k) / max(f_level_psi_pred(:,k)), yi, 'Color', [0.5529    0.6275    0.7961], 'LineWidth', 2);
end

% Customize figure
xlim([0.2, N_L + 1.2])
ylim([0,1.01])
set(gca, 'XTick', 1 : N_L, 'XTickLabel', FACTOR_LABELS)
xlabel('Factor level')
ylabel('Accuracy')
hHline = hline(0.5, 'k--', 'Chance level');
set(gca, 'TickDir', 'out')

hold off

%% Compute quantities of interest
% Summarize the posterior of each approach on the probability scale and on
% the logit scale
for k = 1 : N_L
    fprintf('Posterior median of accuracy for %s is: %g (prob. scale)\n', FACTOR_LABELS{k}, level_psi_prctiles(k,2))
    fprintf('The posterior 95%% CI of accuracy for %s is: [%g, %g] (prob. scale)\n', FACTOR_LABELS{k}, level_psi_prctiles(k,1), level_psi_prctiles(k,3))

     fprintf('Posterior median of accuracy for %s is: %g (logit scale)\n', FACTOR_LABELS{k}, level_alpha_prctiles(k,2))
    fprintf('The posterior 95%% CI of accuracy for %s is: [%g, %g] (logit scale)\n\n', FACTOR_LABELS{k}, level_alpha_prctiles(k,1), level_alpha_prctiles(k,3))
end

% Summarize the posterior of each contrast on the logit scale
for iCntrst = 1 : nCntrsts
    fprintf('Posterior median of contrast %s: %g (logit scale)\n', cntrstLabels{iCntrst}, cntrst_prctiles(iCntrst,2))
    fprintf('The posterior 95%% CI of contrast %s is: [%g, %g] (logit. scale)\n', cntrstLabels{iCntrst}, cntrst_prctiles(iCntrst,1), cntrst_prctiles(iCntrst,3))
    fprintf('Probability that the contrast %s is positive is: %g\n\n', cntrstLabels{iCntrst}, sum(cntrst_vals(:, iCntrst) > 0) / length(cntrst_vals(:, iCntrst)))
end



