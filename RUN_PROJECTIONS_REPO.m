%% RUN_PROJECTIONS_REPO
% Main scenarios for "Reducing the Burden of Cervical Cancer in Zambia in
% the next 40 years: Impact of Enhanced Screening Efforts".
% Run from this repository. Settings below retain the supplied values.
% See README.md for data requirements and preserved model limitations.
clearvars;

%% Settings
pars_local = struct('a2', 0.70, 'a3', 0.05, 'a4', 0.10, ...
    'mu', 0.012, 'Pcum', 0.30);
time_main = 45;
ens = 10;
screenInaccessibleFraction = 0;
calibrationTime = 81;
calibrationEns = 1; % One initial-state simulation, shared by all scenarios.
useParallel = true;

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
outputDir = fullfile(repoRoot, 'results');
outputFile = fullfile(outputDir, 'mainscenarios.mat');

% Columns: [vaccination code, screening code, annual capacity parameter].
% st=0 retains the model's early screening period, then stops screening.
scenarios = [
    0 0 200000   % NOV: no vaccination; screening stops after early period
    1 0 200000   % V: vaccination; screening stops after early period
    1 1 200000   % VT200: vaccination + HPV DNA, capacity parameter 200k
    1 1 400000   % VT400: vaccination + HPV DNA, capacity parameter 400k
    1 2 100000   % VT100: vaccination + VIA, capacity parameter 100k
];
scenarioNames = {'NOV', 'V', 'VT200', 'VT400', 'VT100'};

%% Resolve exact bundled functions and required data
calibFcn = @ABmodelZambia_paper_calibration_only_repo3;
mainFcn = @ABmodelZambia_paper_cal_repo_cases;
assert(exist(func2str(mainFcn), 'file') == 2, 'Missing bundled scenario function.');
assert(isfile(fullfile(repoRoot, 'ageprofile.mat')) || ...
    isfile(fullfile(repoRoot, 'data', 'ageprofile.mat')), ...
    'Place ageprofile.mat in the repository root or data folder.');
fprintf('Using main scenario function: %s\n', func2str(mainFcn));

%% Run exactly one calibration simulation to generate the shared agent pool
assert(exist(func2str(calibFcn), 'file') == 2, ...
    'Missing bundled calibration function.');
assert(isfile(fullfile(repoRoot, 'Pop_1990.xlsx')) || ...
    isfile(fullfile(repoRoot, 'data', 'Pop_1990.xlsx')), ...
    'Place Pop_1990.xlsx in the repository root or data folder.');
fprintf('Generating agenti0 with one calibration simulation...\n');
[~, ~, ~, ~, agenti0, initialCalib] = ...
    calibFcn(pars_local, calibrationTime, calibrationEns);

%% Optional parallel pool; fall back to serial execution if unavailable
runInParallel = false;
if useParallel && license('test', 'Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local');
        end
        runInParallel = true;
    catch ME
        warning('ZambiaScenarios:ParallelUnavailable', ...
            'Could not start the parallel pool; running serially: %s', ME.message);
    end
end

%% Run all five scenarios, retaining incident-case records
nScen = size(scenarios, 1);
emptyResult = struct('name', '', 'sv', [], 'st', [], 'tt', [], ...
    'ensCancer', [], 'ensTested', [], 'SummaryMat', [], ...
    'ensCancer100000', [], 'agenti', [], 'cases', []);
Results = repmat(emptyResult, nScen, 1);
if runInParallel
    parfor k = 1:nScen
        Results(k) = run_one_scenario(mainFcn, pars_local, scenarios(k,:), ...
            scenarioNames{k}, time_main, ens, agenti0, screenInaccessibleFraction);
    end
else
    for k = 1:nScen
        fprintf('Running scenario %d/%d: %s\n', k, nScen, scenarioNames{k});
        Results(k) = run_one_scenario(mainFcn, pars_local, scenarios(k,:), ...
            scenarioNames{k}, time_main, ens, agenti0, screenInaccessibleFraction);
    end
end

%% variables for existing downstream analyses
ensCancerNOV = Results(1).ensCancer;
ensTestedNOV = Results(1).ensTested;
SummaryMatNOV = Results(1).SummaryMat;
ensCancer100000NOV = Results(1).ensCancer100000;
agenti = Results(1).agenti;
casesNOV = Results(1).cases;

ensCancer = Results(2).ensCancer;
ensTested = Results(2).ensTested;
SummaryMat = Results(2).SummaryMat;
ensCancer100000 = Results(2).ensCancer100000;
cases = Results(2).cases;

ensCancerVT200 = Results(3).ensCancer;
ensTestedVT200 = Results(3).ensTested;
SummaryMatVT200 = Results(3).SummaryMat;
ensCancer100000VT200 = Results(3).ensCancer100000;
casesVT200 = Results(3).cases;

ensCancerVT400 = Results(4).ensCancer;
ensTestedVT400 = Results(4).ensTested;
SummaryMatVT400 = Results(4).SummaryMat;
ensCancer100000VT400 = Results(4).ensCancer100000;
casesVT400 = Results(4).cases;

ensCancerVT100 = Results(5).ensCancer;
ensTestedVT100 = Results(5).ensTested;
SummaryMatVT100 = Results(5).SummaryMat;
ensCancer100000VT100 = Results(5).ensCancer100000;
casesVT100 = Results(5).cases;

%% Reproducibility metadata and save
RunMetadata = struct();
RunMetadata.matlabVersion = version;
RunMetadata.calibrationFunction = func2str(calibFcn);
RunMetadata.scenarioFunction = func2str(mainFcn);
RunMetadata.initialStateSource = 'single calibration simulation';
RunMetadata.calibrationWasRun = true;
RunMetadata.calibrationSeed = round(1e6 + 1e3 + calibrationTime);
RunMetadata.runInParallel = runInParallel;
RunMetadata.outputRow = (1:time_main)';
% Requested manuscript labels; review the baseline alignment in README.md.
RunMetadata.requestedYearLabels = (2020:2020+time_main-1)';
RunMetadata.seedByScenario = cell(nScen, 1);
for k = 1:nScen
    RunMetadata.seedByScenario{k} = round(scenarios(k,1)*1e6 + ...
        scenarios(k,2)*1e4 + scenarios(k,3) + (1:ens));
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end
save(outputFile, ...
    'Results', 'RunMetadata', ...
    'pars_local', 'time_main', 'ens', 'scenarios', 'scenarioNames', ...
    'screenInaccessibleFraction', 'calibrationTime', 'calibrationEns', ...
    'initialCalib', ...
    'agenti0', 'agenti', ...
    'ensCancerNOV', 'ensTestedNOV', 'SummaryMatNOV', 'ensCancer100000NOV', ...
    'ensCancer', 'ensTested', 'SummaryMat', 'ensCancer100000', ...
    'ensCancerVT200', 'ensTestedVT200', 'SummaryMatVT200', 'ensCancer100000VT200', ...
    'ensCancerVT400', 'ensTestedVT400', 'SummaryMatVT400', 'ensCancer100000VT400', ...
    'ensCancerVT100', 'ensTestedVT100', 'SummaryMatVT100', 'ensCancer100000VT100', ...
    'casesNOV', 'cases', 'casesVT200', 'casesVT400', 'casesVT100', '-v7.3');
fprintf('Saved results to %s\n', outputFile);

%% Local helper
function result = run_one_scenario(mainFcn, pars, scenario, name, ...
        timeSteps, ens, agenti0, inaccessibleFraction)
    result = struct('name', name, 'sv', scenario(1), 'st', scenario(2), ...
        'tt', scenario(3), 'ensCancer', [], 'ensTested', [], 'SummaryMat', [], ...
        'ensCancer100000', [], 'agenti', [], 'cases', []);
    [result.ensCancer, result.ensTested, result.SummaryMat, ...
        result.ensCancer100000, result.agenti, result.cases] = ...
        mainFcn(pars, scenario(1), scenario(2), scenario(3), ...
            timeSteps, ens, false, agenti0, inaccessibleFraction);
end
