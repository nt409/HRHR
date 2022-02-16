# HRHR project

## Model

Model based on those of Hobbelen and Elderfield. Two high risk fungicides - we
model how the pathogen population becomes more resistant as the fungicides are
applied over a number of years.

## Code

The model is implemented in python.

In the work 'Optimal resistance management for mixtures of high-risk
fungicides: robustness to the initial frequency of resistance and pathogen
sexual reproduction', we describe the model found in `src\model`.

The most important file in this folder is `simulator.py`, which contains the
classes `RunSingleTactic` and `RunGrid` which are the main ways we test tactic
performance. These have docstrings which describe their use and how to
configure them.

### Scans

We ran three scans over different parameter values/scenarios in the work:

- Robustness to variation in parameter values; see `src\param_scan`
- Effect of between-season sexual reproduction; see `src\sr_scan`
- Mixtures vs Alternations (with between-season sexual reproduction); see
  `src\alternation_scan`

### Plotting

Some plotting functions are found in `src\plotting`, although you may prefer to
write your own custom plotting functions. Most of these are called in various
scripts in `src\figs` or in the notebooks.
