# Bayesian BCI performance

In this repository you can find the MATLAB code that accompanies the paper "Beyond p-values in the evaluation of brain-computer interfaces" (preprint coming soon).
The implemented models can be used to analyze experimentally obtained brain-computer interface (BCI) classification accuracies.
The code can be used to substitute the classical hypothesis testing with Bayesian estimation of hierarchical models.

## Requirements

The code relies on the freely available WinBUGS software to perform Bayesian inference using Markov chain Monte Carlo (MCMC) sampling.
The WinBUGS software is distributed separately and can be obtained at:
http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/

It should also be possible to use the implemented models with other MCMC software based on the BUGS language, such as OpenBUGS or JAGS,
but this has not yet been tested. Whereas OpenBUGS can be interfaced from MATLAB using the same interface which we used for WinBUGS
(i.e. MATBUGS), to interface with JAGS it is necessary to use the MATJAGS interface available at:
http://psiexp.ss.uci.edu/research/programs_data/jags/

## Usage

The repository contains four models: one example model for simple normal data and three hierarchical models for different BCI experimental designs.

The models for BCI performance are:
**Model 1**
  * appropriate for an experiment where a group of subjects has used a BCI
  * for each subject the number of correct trials and total number of trials is needed
  * this data would typically be analyzed with a single sample t-test
  
**Model 2**
  * appropriate for an experiment where we wish to estimate an association between a subject-specific continuous variable (could be observed, like age, or experimentally set, like amount of time that user was given to learn to use the BCI)
  * the data needs to contain the values of the subject-specific explanatory variable as well as the number of correct trials and the total number of trials per subject 
  * this data would typically be analyzed with a linear regression

**Model 3**
  - appropriate for an experiment where each subject has tried multiple BCI approaches (could be different experimental approaches, like different stimuli, or different computational approaches, like different classifiers)
  - the data consists of observations defined by a subject code, condition code, and the number of correct trials and the total number of trials for the subject-condition combination
  - this data would typically be analyzed with a matched t-test (if there are two conditions), or a repeated measures ANOVA (for multiple conditions)
  

The folder of each model has the following contents:
./data/
simulate_modelN_data.m
winbugs-modelN.txt
run_modelN.m
inspect_modelN_results.m

with N = 0, 1, 2, 3.

The data folder contains the example datasets for the model -- one simulated dataset with known true parameter values
and one real dataset (as described in the reference paper). The .csv files are used to input the data, and the example files show
the required formatting of the data.

The **simulate_modelN_data** script is used to generate synthetic datasets, which can in turn be used to test the inference procedure.
This script does not need to be used to analyze real datasets.

The **winbugs-modelN** text file is the description of the statistical model in the BUGS modeling language.

The **run_modelN** script is used to perform Bayesian inference via MCMC sampling in WinBUGS software.
Running this script will open the WinBUGS program and the results of the inference can also be inspected there.
To finish the analysis, WinBUGS needs to be closed, after which MATLAB will continue executing the script.

The **inspect_modelN_results** script will take the results generated by *run_modelN* script and report the
results of the inference, i.e. different properties of the posterior distribution of the model's parameters.

Before running a script, make sure that the options (such as path variables) at the top of the script are set correctly.
For example, the variable PATH_BUGS in **run_modelN** needs to be set to the path of the WinBUGS installation.


## Troubleshooting

The MATLAB-WinBUGS interface is implemented by the matbugs.m function in /external/matbugs folder.
However, this function did not work when we tested it on Windows 7 with MATLAB R2011a (due to an issue with a call to the dos() function).
If that happens to be the case with your configuration, try to replace the function with the matbugs_alternative.m function in the same folder.

The code should also be able to run directly with OpenBUGS but we have not yet tested this.

## Credits

The code is maintained by Filip Melinscak (filip.melinscak@bitbrain.es).

## License

The code is released under the GNU GPL-3 license.