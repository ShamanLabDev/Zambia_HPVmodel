 Zambia cervical cancer model

## Setup

Use MATLAB with Statistics and Machine Learning Toolbox. Parallel Computing Toolbox is optional.

Place these files in `data/`:

- `Pop_1990.xlsx` — required for calibration and projections.
- `ageprofile.mat` — required for projections.

Set MATLAB's current folder to this repository.

## Calibration

Edit the parameter grids and ensemble count in `RUN_CALIBRATION_GRID_REPO.m`, then run:

RUN_CALIBRATION_GRID_REPO


Calibration results saved in `results/`.

## Projections

Set `pars_local` and the scenario settings in `RUN_PROJECTIONS_REPO.m`, then run:

RUN_PROJECTIONS_REPO

The script runs one calibration simulation to generate `agenti0`, then uses that pool for all five projection scenarios. Results are saved to `results/mainscenarios.mat`.

To use parameters from the calibration grid, copy `pars_best` into `pars_local` before running projections. The grid does not need to be rerun when using an existing parameter set.
