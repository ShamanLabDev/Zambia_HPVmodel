function [ensCancer, ensTested, SummaryMat, ensCancer100000, agentstable,cases] = ...
    ABmodelZambia_paper_cal_repo_cases(pars, sv, st, tt, timeSteps, ens, doPlot, agenti0, screenInaccessibleFraction)
% ABMODELZAMBIA_PAPER_CAL_REPO_CASES
% Repository-ready version of the Zambia HPV/cervical cancer IBM.
%%
% Required external files/functions (expected on MATLAB path or in repo):
%   - ageprofile.mat        (contains age19F, age19M)
%
% Inputs
%   pars      struct with fields: mu, Pcum, a2, a3, a4
%   sv        vaccination scenario code
%   st        screening scenario code (0 = no future testing, 1 = HPV DNA, 2 = VIA)
%   tt        annual test capacity parameter
%   timeSteps number of yearly simulation steps
%   ens       number of ensemble simulations
%   doPlot    logical, optional
%   agenti0   table with initialized agent states from Calibration
%   screenInaccessibleFraction  persistent fraction of the population that
%              can never be selected for screening; default = 0
%
% Outputs
%   ensCancer         yearly cervical cancer counts mapped to real population
%   ensTested         yearly number tested mapped to real population
%   SummaryMat        summary table from the last ensemble run
%   ensCancer100000   yearly cancer incidence per 100,000 women
%   agentstable       final agent table from the last ensemble run
%   cases             cell array (one cell per ensemble); each row is
%                     [age at cancer onset, years since HPV acquisition]
% Screening-access sensitivity:
%   Each individual receives a persistent Uniform(0,1) ScreenAccessScore at
%   initialization or birth. Screening is limited to women with
%   ScreenAccessScore <= 1-screenInaccessibleFraction. Thus, a fixed subgroup
%   remains outside the screening-reachable population throughout life.
%
%% -------------------------------------------------------------------------
    narginchk(7, 9);
    if nargin == 7 && istable(doPlot)
        agenti0 = doPlot;
        doPlot = false;
    elseif nargin < 8
        error(['ABmodelZambia_paper_cal_repo requires either 7 inputs ' ...
               '(with agenti0 passed in place of doPlot), 8, or 9 inputs.']);
    end
    if nargin < 7 || isempty(doPlot)
        doPlot = false;
    end
    if nargin < 9 || isempty(screenInaccessibleFraction)
        screenInaccessibleFraction = 0;
    end
    validateattributes(screenInaccessibleFraction, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1}, ...
        mfilename, 'screenInaccessibleFraction', 9);
    screenReachableFraction = 1 - screenInaccessibleFraction;
    validateattributes(pars, {'struct'}, {'scalar'}, mfilename, 'pars', 1);
    requiredPars = {'mu','Pcum','a2','a3','a4'};
    for k = 1:numel(requiredPars)
        if ~isfield(pars, requiredPars{k})
            error('Missing pars.%s', requiredPars{k});
        end
    end
    validateattributes(sv, {'numeric'}, {'real','finite','scalar','integer','>=',0,'<=',18});
    validateattributes(st, {'numeric'}, {'real','finite','scalar','integer','>=',0,'<=',2});
    validateattributes(tt, {'numeric'}, {'real','finite','scalar','nonnegative'});
    validateattributes(timeSteps, {'numeric'}, {'real','finite','scalar','integer','>=',6});
    validateattributes(ens, {'numeric'}, {'real','finite','scalar','integer','>=',1});
    validateattributes(doPlot, {'logical','numeric'}, {'scalar'});
    if ~istable(agenti0)
        error('agenti0 must be a MATLAB table.');
    end
    for k = 1:numel(requiredPars)
        validateattributes(pars.(requiredPars{k}), {'numeric'}, ...
            {'real','finite','scalar','nonnegative'}, mfilename, requiredPars{k});
    end
    if pars.mu > 1 || pars.Pcum > 1 || ...
            any([pars.a2, pars.a2*pars.a3, pars.a2*pars.a3*pars.a4] > 1)
        error('ZambiaScenarios:InvalidProbability', ...
            'mu, Pcum, a2, a2*a3, and a2*a3*a4 must be in [0, 1].');
    end
    requiredAgentFields = {'Age','Gender','HIVStatus','VaccineStatus', ...
        'InfectionStatus','TestedStatus','LesionStatus','CancerStatus', ...
        'YearsSinceTested','IsAlive','YearsSinceHPV','YearsSinceLesions', ...
        'YearsSinceCancer'};
    if ~all(ismember(requiredAgentFields, agenti0.Properties.VariableNames))
        error('ZambiaScenarios:InvalidInitialState', ...
            'agenti0 is missing one or more required agent-state columns.');
    end
    callerRng = rng;
    restoreRng = onCleanup(@() rng(callerRng)); 
    repoRoot = fileparts(mfilename('fullpath'));
    [age19F, age19M] = load_age_profiles(repoRoot);
    agenti = agenti0;
    if ismember('IsAlive', agenti.Properties.VariableNames)
        agenti = agenti(agenti.IsAlive == 1, :);
    end
    populationSize = 50000;   % total agents (male + female)
    realSizeM      = 8994252; % males in 2019 projections
    realSizeF      = 8984359; % females in 2019 projections
    realSizeTOT    = realSizeF + realSizeM;
    initial_growth_rate = 0.028;
    decrease_rate       = 1/100;
    populationSizeY     = growth_rate_paper(populationSize, initial_growth_rate, decrease_rate, timeSteps);
    ageGroups = [14, 24, 34, 44, 54];
    agesmoothRF = smoothdata(repelem(age19F/5, 5), 'movmean', 5);
    agesmoothRM = smoothdata(repelem(age19M/5, 5), 'movmean', 5);
    normaliz    = populationSize / realSizeTOT;
    agesmoothF  = round(agesmoothRF * normaliz);
    agesmoothM  = round(agesmoothRM * normaliz);
    agesmooth   = agesmoothF + agesmoothM;
    if sum(agesmooth) < populationSize
        error('ZambiaScenarios:InsufficientAgeProfile', ...
            'The rounded age profile supplies %d agents; %d are required.', ...
            sum(agesmooth), populationSize);
    end
    requiredAges = find(agesmooth > 0) - 1;
    if ~all(ismember(requiredAges, agenti.Age))
        error('ZambiaScenarios:MissingAgeDonor', ...
            'The living initial-state table lacks donors for a required exact age.');
    end
    deathProbability = storedata();
    wASR84 = [12000 10000 9000 9000 8000 8000 6000 6000 6000 ...
              6000 5000 4000 4000 3000 2000 1000 500];
    probabilityLesion = pars.mu;
    Pcum              = pars.Pcum;
    a2                = pars.a2;
    a3                = pars.a3;
    a4                = pars.a4;
    s2   = a2;
    s3   = a2 * a3;
    s4   = a2 * a3 * a4;
    SHPV = [0, s2, s3, s4, 0, 0];
    incHIVM = [0, 0.0006, 0.0006, 0.0006, 0, 0];
    incHIVF = [0, 0.0063, 0.0063, 0.0063, 0, 0];
    HIVincFunction = @(Age, Gender) ...
        (Gender == 0) .* ( ...
        (Age <= ageGroups(1))                     * incHIVM(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * incHIVM(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * incHIVM(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * incHIVM(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * incHIVM(5) + ...
        (Age > ageGroups(5))                       * incHIVM(6)) + ...
        (Gender == 1) .* ( ...
        (Age <= ageGroups(1))                     * incHIVF(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * incHIVF(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * incHIVF(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * incHIVF(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * incHIVF(5) + ...
        (Age > ageGroups(5))                       * incHIVF(6));
    RR_acq    = 1.7;
    RR_lesion = 2.4;
    HR_clear  = 0.7;
    mHIV      = 0.01;
    prevalenceValues0F = [0, 0.35, 0.3, 0.2, 0.15, 0.15];
    prevalenceValues0M = [0, 0.3, 0.4, 0.35, 0.25, 0.25];
    BW = [0, 1.5, 1, 0.8, 0.5, 0.25];
    BM = [0, 1.7, 2, 2, 1.3, 1];
    PartnerFunction = @(Age, Gender) ...
        (Gender == 0) .* ( ...
        (Age <= ageGroups(1))                     * BM(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * BM(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * BM(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * BM(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * BM(5) + ...
        (Age > ageGroups(5))                       * BM(6)) + ...
        (Gender == 1) .* ( ...
        (Age <= ageGroups(1))                     * BW(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * BW(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * BW(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * BW(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * BW(5) + ...
        (Age > ageGroups(5))                       * BW(6));
    susceptibility = @(Age) ...
        (Age <= ageGroups(1))                     * SHPV(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * SHPV(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * SHPV(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * SHPV(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * SHPV(5) + ...
        (Age > ageGroups(5))                       * SHPV(6);
    recoveryrate1 = @(YearsSinceHPV) ...
        (YearsSinceHPV == 1) * 0.6 + ...
        (YearsSinceHPV == 2) * 0.5 + ...
        (YearsSinceHPV == 3) * 0.2 + ...
        (YearsSinceHPV == 4) * 0.2;
    Py = 1 - (1 - Pcum)^(1/15);
    CancerRatefromLesion = @(YearsSinceLesions) (YearsSinceLesions > 5) * Py;
    DeathCancer          = @(YearsSinceCancer) ...
        (YearsSinceCancer <= 5) * 0.12 + ...
        (YearsSinceCancer <= 10 & YearsSinceCancer > 5) * 0.05;
    sensitivity       = 0.98;
    sensitivityVIA    = 0.35;
    efficacyHIV       = 0.63;
    efficacyNEG       = 0.85;
    ratiotreated      = 0.9;
    ratiotreatedVIA   = 0.73;
    probabilityTreated    = sensitivity    * ratiotreated;
    probabilityTreatedVIA = sensitivityVIA * ratiotreatedVIA;
    Ppositive_after_treatment = 0.1;
    multi = realSizeF / sum(agesmoothF);
    hpv_last10_all    = zeros(ens, 6);
    lesion_last10_all = zeros(ens, 1);
    cancer_last10_all = zeros(ens, 1);
    ensCancer       = NaN(timeSteps, ens);
    ensTested       = NaN(timeSteps, ens);
    ensCancer100000 = NaN(timeSteps, ens);
    SummaryMat      = table();
    agentstable     = table();
    % Incident cancer records for each stochastic ensemble.
    % Each row of cases{nsim} is:
    %   [age at cancer onset, years since HPV acquisition]
    cases = cell(ens,1);
    for nsim = 1:ens
        baseSeed = sv * 1e6 + st * 1e4 + tt + nsim;
        rng(round(baseSeed), 'twister');
        agents = struct( ...
            'Age', 0, ...
            'Gender', 0, ...
            'HIVStatus', 0, ...
            'VaccineStatus', 0, ...
            'InfectionStatus', 0, ...
            'TestedStatus', 0, ...
            'LesionStatus', 0, ...
            'CancerStatus', 0, ...
            'YearsSinceTested', 100, ...
            'IsAlive', 1, ...
            'YearsSinceHPV', 0, ...
            'YearsSinceLesions', 0, ...
            'YearsSinceCancer', 0, ...
            'TreatmentStatus', 0, ...
            'ScreenAccessScore', 0);
        ages = 0:84;
        age_vector = repelem(ages, agesmooth);
        age_vector = age_vector(randperm(length(age_vector)));
        [agents, agentstable0] = initialization(populationSize, age_vector, agenti, agents);
        if nsim == 1
            femaleInitial = agentstable0.Gender == 1 & agentstable0.IsAlive == 1;
            realizedReach = mean(agentstable0.ScreenAccessScore(femaleInitial) <= screenReachableFraction);
            fprintf('Screening sensitivity: %.1f%% permanently inaccessible; realized reachable women = %.1f%%\n', ...
                100 * screenInaccessibleFraction, 100 * realizedReach);
        end
        Pepfar = 0.58 * [130000, 195000, 247000, 268300, 266700];
        yearlytesttargetNEG = zeros(timeSteps, 1);
        yearlytesttargetHIV = zeros(timeSteps, 1);
        leftovers      = max(0, tt - Pepfar(5));
        rampupnumberHIV = round(Pepfar / multi);
        rampupnumberNEG = round(leftovers / multi);
        if st == 1 || st == 2
            yearlytesttargetHIV(1:5) = round(rampupnumberHIV * 0.6);
            yearlytesttargetNEG(1:5) = round(rampupnumberHIV * 0.4);
            yearlytesttargetNEG(6)   = yearlytesttargetNEG(5) + rampupnumberNEG;
            yearlytesttargetHIV(6)   = yearlytesttargetHIV(5);
        elseif st == 0
            yearlytesttargetHIV(1:5) = round(rampupnumberHIV * 0.6);
            yearlytesttargetNEG(1:5) = round(rampupnumberHIV * 0.4);
            yearlytesttargetNEG(6:timeSteps) = 0;
            yearlytesttargetHIV(6:timeSteps) = 0;
        end
        yearlytesttargetTOT = yearlytesttargetHIV + yearlytesttargetNEG;
        yearlytesttargetTOT(1:6) = yearlytesttargetHIV(1:6) + yearlytesttargetNEG(1:6);
        for iYear = 7:timeSteps
            yearlytesttargetHIV(iYear) = round(yearlytesttargetHIV(iYear-1) + 0.05 * yearlytesttargetHIV(iYear-1));
            yearlytesttargetNEG(iYear) = round(yearlytesttargetNEG(iYear-1) + 0.05 * yearlytesttargetNEG(iYear-1));
            yearlytesttargetTOT(iYear) = yearlytesttargetNEG(iYear) + yearlytesttargetHIV(iYear);
        end
        for i = populationSizeY(1)+1:populationSizeY(timeSteps)
            agents(i).Age = NaN;
            agents(i).Gender = NaN;
            agents(i).HIVStatus = NaN;
            agents(i).TestedStatus = NaN;
            agents(i).LesionStatus = NaN;
            agents(i).InfectionStatus = NaN;
            agents(i).CancerStatus = NaN;
            agents(i).VaccineStatus = NaN;
            agents(i).YearsSinceTested = NaN;
            agents(i).YearsSinceLesions = NaN;
            agents(i).YearsSinceHPV = NaN;
            agents(i).IsAlive = NaN;
            agents(i).YearsSinceCancer = NaN;
            agents(i).TreatmentStatus = NaN;
            agents(i).ScreenAccessScore = NaN;
        end
        mat          = zeros(timeSteps, 24);
        hpv_by_year  = nan(timeSteps, 6);
        ASR84_by_year = NaN(timeSteps, 1);
        Ncured       = 0;
        female0 = agentstable0(agentstable0.Gender == 1 & agentstable0.IsAlive == 1, :);
        hpv_by_year(1,1) = mean(female0.InfectionStatus(female0.Age <= 14));
        hpv_by_year(1,2) = mean(female0.InfectionStatus(female0.Age >= 15 & female0.Age <= 24));
        hpv_by_year(1,3) = mean(female0.InfectionStatus(female0.Age >= 25 & female0.Age <= 34));
        hpv_by_year(1,4) = mean(female0.InfectionStatus(female0.Age >= 35 & female0.Age <= 44));
        hpv_by_year(1,5) = mean(female0.InfectionStatus(female0.Age >= 45 & female0.Age <= 54));
        hpv_by_year(1,6) = mean(female0.InfectionStatus(female0.Age >= 55));
        hpv_by_year(1, isnan(hpv_by_year(1,:))) = 0;
        Tcured = agentstable0(1,:);
        Tcancer = table();
        cancerDeathCounter = 0;
        elig = zeros(timeSteps, 5);
        elig(1,:) = 0;
        cases_nsim = zeros(0, 2);
        for t = 2:timeSteps
            incCancerAge84 = zeros(17,1);
            incCancer      = 0;
            incCancerHIV   = 0;
            incCancerNOHIV = 0;
            incCancer30    = 0;
            incCancer40    = 0;
            incCancer50    = 0;
            incCancer60    = 0;
            incCancerDeath = 0;
            agentstable = struct2table(agents);
            state0      = agents;
            nsexvec     = poissrnd(PartnerFunction(agentstable.Age, agentstable.Gender));
            cover = 0.9 * [0.99, 0.73, 0.44, 0.66, 0.56];
            if t < 7 && sv > 0
                eligible_agents_vac = find(agentstable.Gender == 1 & agentstable.Age == 13);
                vyear = round(cover(t-1) * size(eligible_agents_vac, 1));
                if ~isempty(eligible_agents_vac) && vyear > 0
                    selected_idV = randsample(eligible_agents_vac, min(vyear, numel(eligible_agents_vac)));
                    [agents(selected_idV).VaccineStatus] = deal(1);
                end
            elseif t >= 7 && sv > 0
                switch sv
                    case 1
                        agents = vaccinate_agents(0.9 * 0.73, agentstable.Gender == 1, agentstable, agents);
                    case 2
                        agents = vaccinate_agents(0.9 * 0.91, agentstable.Gender == 1, agentstable, agents);
                    case 3
                        agents = vaccinate_agents(0.9 * 0.46, agentstable.Gender < 2, agentstable, agents);
                    case 4
                        agents = vaccinate_agents(0.9 * 0.91, agentstable.Gender < 2, agentstable, agents);
                    case 5
                        agents = vaccinate_agents(0.9 * 0.73, agentstable.Gender < 2, agentstable, agents);
                    case 6
                        agents = vaccinate_agents(0.9 * 0.46, agentstable.Gender == 1, agentstable, agents);
                    case 7
                        agents = vaccinate_agents(0.9 * 0.64, agentstable.Gender == 1, agentstable, agents);
                    case 8
                        agents = vaccinate_agents(0.9 * 0.64, agentstable.Gender < 2, agentstable, agents);
                    case 9
                        agents = vaccinate_agents(0.9 * 0.55, agentstable.Gender == 1, agentstable, agents);
                    case 10
                        agents = vaccinate_agents(0.9 * 0.55, agentstable.Gender < 2, agentstable, agents);
                    case 11
                        agents = vaccinate_agents(0.9 * 0.82, agentstable.Gender == 1, agentstable, agents);
                    case 12
                        agents = vaccinate_agents(0.9 * 0.82, agentstable.Gender < 2, agentstable, agents);
                    case 13
                        agents = vaccinate_agents(0.9 * 0.37, agentstable.Gender == 1, agentstable, agents);
                    case 14
                        agents = vaccinate_agents(0.9 * 0.37, agentstable.Gender < 2, agentstable, agents);
                    case 15
                        agents = vaccinate_agents(0.9 * 0.28, agentstable.Gender == 1, agentstable, agents);
                    case 16
                        agents = vaccinate_agents(0.9 * 0.28, agentstable.Gender < 2, agentstable, agents);
                    case 17
                        agents = vaccinate_agents(0.9 * 0.19, agentstable.Gender == 1, agentstable, agents);
                    case 18
                        agents = vaccinate_agents(0.9 * 0.19, agentstable.Gender < 2, agentstable, agents);
                end
            end
            for i = 1:populationSizeY(t)
                agents(i).Age = agents(i).Age + 1;
                if agents(i).Age == 85 || agents(i).IsAlive == 0 || i > populationSizeY(t-1)
                    agents(i).Age = 0;
                    agents(i).Gender = double(rand() < 0.5);
                    agents(i).HIVStatus = double(rand() < 0.005);
                    agents(i).InfectionStatus = 0;
                    agents(i).LesionStatus = 0;
                    agents(i).VaccineStatus = 0;
                    agents(i).YearsSinceLesions = 0;
                    agents(i).YearsSinceHPV = 0;
                    agents(i).CancerStatus = 0;
                    agents(i).TestedStatus = 0;
                    agents(i).YearsSinceTested = 100;
                    agents(i).IsAlive = 1;
                    agents(i).YearsSinceCancer = 0;
                    agents(i).TreatmentStatus = 0;
                    % Persistent lifetime screening-access score. Women with
                    % scores above screenReachableFraction are never eligible.
                    agents(i).ScreenAccessScore = rand();
                    continue;
                end
                age = agents(i).Age;
                deathP = deathProbability(max(round(age),1)) + mHIV * agents(i).HIVStatus;
                if rand() < deathP
                    agents(i).IsAlive = 0;
                    continue;
                end
                if agents(i).Age >= 15
                    if rand() < HIVincFunction(agents(i).Age, agents(i).Gender) + agents(i).HIVStatus
                        agents(i).HIVStatus = 1;
                    end
                    if agents(i).Gender == 0 && agents(i).VaccineStatus == 0
                        if agents(i).InfectionStatus == 0
                            nsex = nsexvec(i);
                            if nsex > 0
                                eligible_agents_idx = find(agentstable.Gender == 1 & (age - agentstable.Age) <= 5 & ...
                                    (age - agentstable.Age) >= -1 & agentstable.Age >= 15);
                                if ~isempty(eligible_agents_idx)
                                    selected_idx = eligible_agents_idx(randi(numel(eligible_agents_idx), nsex, 1));
                                    m = sum([state0(selected_idx).InfectionStatus]);
                                    beta = susceptibility(age);
                                    p_any0 = 1 - (1 - beta)^m;
                                    p_any0 = min(max(p_any0, 0), 1);
                                    if agents(i).HIVStatus == 1
                                        p_any = 1 - (1 - p_any0)^RR_acq;
                                    else
                                        p_any = p_any0;
                                    end
                                    if rand() < p_any
                                        agents(i).InfectionStatus = 1;
                                        agents(i).YearsSinceHPV   = 1;
                                    end
                                end
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 0
                            pippo = rand();
                            p1_0 = probabilityLesion;
                            if agents(i).HIVStatus == 1
                                p1 = 1 - (1 - p1_0)^RR_lesion;
                            else
                                p1 = p1_0;
                            end
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if pippo < p1
                                agents(i).LesionStatus = 1;
                                agents(i).YearsSinceLesions = 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            elseif pippo < p1 + pclr
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1
                            pippo = rand();
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if pippo < pclr
                                agents(i).LesionStatus = 0;
                                agents(i).YearsSinceLesions = 0;
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        end
                    elseif agents(i).Gender == 1
                        if agents(i).InfectionStatus == 0 && agents(i).VaccineStatus == 0
                            nsex = nsexvec(i);
                            if nsex > 0
                                eligible_agents_idx = find(agentstable.Gender == 0 & (agentstable.Age - age) <= 5 & ...
                                    (agentstable.Age - age) >= -1 & agentstable.Age >= 15);
                                if ~isempty(eligible_agents_idx)
                                    selected_idx = eligible_agents_idx(randi(numel(eligible_agents_idx), nsex, 1));
                                    m = sum([state0(selected_idx).InfectionStatus]);
                                    beta = susceptibility(age);
                                    p_any0 = 1 - (1 - beta)^m;
                                    p_any0 = min(max(p_any0,0),1);
                                    if agents(i).HIVStatus == 1
                                        p_any = 1 - (1 - p_any0)^RR_acq;
                                    else
                                        p_any = p_any0;
                                    end
                                    if rand() < p_any
                                        agents(i).InfectionStatus = 1;
                                        agents(i).YearsSinceHPV = 1;
                                    end
                                end
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 0
                            pippo = rand();
                            p1_0 = probabilityLesion;
                            if agents(i).HIVStatus == 1
                                p1 = 1 - (1 - p1_0)^RR_lesion;
                            else
                                p1 = p1_0;
                            end
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if pippo < p1
                                agents(i).LesionStatus = 1;
                                agents(i).YearsSinceLesions = 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            elseif pippo < p1 + pclr
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1 && agents(i).CancerStatus == 0
                            pippo = rand();
                            p2 = CancerRatefromLesion(agents(i).YearsSinceLesions);
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if pippo < p2
                                agents(i).CancerStatus = 1;
                                incCancer = incCancer + 1;
                                incCancerHIV = incCancerHIV + agents(i).HIVStatus;
                                incCancerNOHIV = incCancerNOHIV + (1 - agents(i).HIVStatus);
                                band = ageBand84(agents(i).Age);
                                % Record incident cancer at the exact
                                % transition from CancerStatus 0 -> 1.
                                % Col 1: age at cancer onset
                                % Col 2: years since HPV acquisition
                                cases_nsim(end+1,:) = [agents(i).Age, agents(i).YearsSinceHPV]; %#ok<AGROW>
                                if ~isnan(band)
                                    incCancerAge84(band) = incCancerAge84(band) + 1;
                                end
                                incCancer30 = incCancer30 + (agents(i).Age >= 25 && agents(i).Age < 35);
                                incCancer40 = incCancer40 + (agents(i).Age >= 35 && agents(i).Age < 45);
                                incCancer50 = incCancer50 + (agents(i).Age >= 45 && agents(i).Age < 55);
                                incCancer60 = incCancer60 + (agents(i).Age >= 55);
                                agents(i).YearsSinceCancer = 1;
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            elseif pippo < p2 + pclr
                                agents(i).LesionStatus = 0;
                                agents(i).YearsSinceLesions = 0;
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1 && agents(i).CancerStatus == 1
                            pippo = rand();
                            p3 = DeathCancer(agents(i).YearsSinceCancer + mHIV * agents(i).HIVStatus);
                            if pippo < p3
                                agents(i).IsAlive = 0;
                                cancerDeathCounter = cancerDeathCounter + 1;
                                Tcancer(cancerDeathCounter, :) = struct2table(agents(i));
                            else
                                agents(i).YearsSinceCancer = agents(i).YearsSinceCancer + 1;
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                            end
                        end
                        agents(i).YearsSinceTested = agents(i).YearsSinceTested + 1;
                    end
                end
            end
            agentstable = struct2table(agents);
            eligibleagentsHIVneg = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 5 & ...
                agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & ...
                agentstable.HIVStatus == 0 & ...
                agentstable.ScreenAccessScore <= screenReachableFraction);
            eligibleagentsHIVpos = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 3 & ...
                agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & ...
                agentstable.HIVStatus == 1 & ...
                agentstable.ScreenAccessScore <= screenReachableFraction);
            needed = length(eligibleagentsHIVpos) + length(eligibleagentsHIVneg);
            TestAnnoHIV = min(yearlytesttargetHIV(t), length(eligibleagentsHIVpos));
            left = max(0, yearlytesttargetHIV(t) - length(eligibleagentsHIVpos));
            yearlytesttargetNEG(t) = yearlytesttargetNEG(t) + left;
            yearlytesttargetTOT(t) = yearlytesttargetHIV(t) + yearlytesttargetNEG(t);
            if TestAnnoHIV > 0
                testwinnersID_HIV = randperm(length(eligibleagentsHIVpos), TestAnnoHIV);
            else
                testwinnersID_HIV = [];
            end
            for iWin = testwinnersID_HIV
                yy = eligibleagentsHIVpos(iWin);
                agents(yy).TestedStatus = agents(yy).TestedStatus + 1;
                agents(yy).YearsSinceTested = 0;
                if t < 6 || (t > 5 && st == 2 && agents(yy).LesionStatus == 1 && agents(yy).CancerStatus == 0)
                    if rand() < probabilityTreatedVIA * efficacyHIV
                        agents(yy).TreatmentStatus = 1;
                        agents(yy).LesionStatus = 0;
                        agents(yy).YearsSinceLesions = 0;
                        Ncured = Ncured + 1;
                        Tcured = vertcat(Tcured, agentstable(yy,:));
                        if rand() < (1 - Ppositive_after_treatment)
                            agents(yy).InfectionStatus = 0;
                        end
                    end
                elseif t > 5 && st == 1 && agents(yy).CancerStatus == 0
                    if rand() < probabilityTreated * efficacyHIV
                        agents(yy).TreatmentStatus = 1;
                        agents(yy).LesionStatus = 0;
                        agents(yy).YearsSinceLesions = 0;
                        Ncured = Ncured + 1;
                        Tcured = vertcat(Tcured, agentstable(yy,:));
                        if rand() < (1 - Ppositive_after_treatment)
                            agents(yy).InfectionStatus = 0;
                        end
                    end
                end
            end
            agentstable = struct2table(agents);
            eligibleagentsHIVgen = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 5 & ...
                agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & ...
                agentstable.ScreenAccessScore <= screenReachableFraction);
            TestAnnoNEG = min(yearlytesttargetNEG(t), length(eligibleagentsHIVgen));
            if TestAnnoNEG > 0
                testwinnersID_NEG = randperm(length(eligibleagentsHIVgen), TestAnnoNEG);
            else
                testwinnersID_NEG = [];
            end
            for iWin = testwinnersID_NEG
                idx = eligibleagentsHIVgen(iWin);
                agents(idx).TestedStatus = agents(idx).TestedStatus + 1;
                agents(idx).YearsSinceTested = 0;
                if t < 6 || (st == 2 && agents(idx).LesionStatus == 1 && agents(idx).CancerStatus == 0)
                    if rand() < probabilityTreatedVIA * agents(idx).HIVStatus * efficacyHIV + ...
                            probabilityTreatedVIA * (1 - agents(idx).HIVStatus) * efficacyNEG
                        agents(idx).TreatmentStatus = 1;
                        agents(idx).LesionStatus = 0;
                        agents(idx).YearsSinceLesions = 0;
                        Ncured = Ncured + 1;
                        Tcured = vertcat(Tcured, agentstable(idx,:));
                        if rand() < (1 - Ppositive_after_treatment)
                            agents(idx).InfectionStatus = 0;
                        end
                    end
                elseif t > 5 && st == 1 && agents(idx).CancerStatus == 0
                    agents(idx).TreatmentStatus = 1;
                    if rand() < probabilityTreated * agents(idx).HIVStatus * efficacyHIV + ...
                            probabilityTreated * (1 - agents(idx).HIVStatus) * efficacyNEG
                        agents(idx).TreatmentStatus = 1;
                        agents(idx).LesionStatus = 0;
                        agents(idx).YearsSinceLesions = 0;
                        Ncured = Ncured + 1;
                        Tcured = vertcat(Tcured, agentstable(idx,:));
                        if rand() < (1 - Ppositive_after_treatment)
                            agents(idx).InfectionStatus = 0;
                        end
                    end
                end
            end
            agentstable = struct2table(agents);
            atnormHIVneg = sum(agentstable.Gender == 1 & agentstable.YearsSinceTested <= 5 & agentstable.Age >= 25 & ...
                agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 0);
            atnormHIVpos = sum(agentstable.Gender == 1 & agentstable.YearsSinceTested <= 3 & agentstable.Age >= 25 & ...
                agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 1);
            tottargetHIVneg = sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & ...
                agentstable.IsAlive == 1 & agentstable.HIVStatus == 0);
            tottargetHIVpos = sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & ...
                agentstable.IsAlive == 1 & agentstable.HIVStatus == 1);
            sumValues = zeros(1, 24);
            agents2 = agents([agents.IsAlive] == 1 & [agents.Gender] == 1);
            denom = sum(agentstable.Gender == 1);
            for iAlive = 1:length(agents2)
                sumValues(1) = sumValues(1) + agents2(iAlive).InfectionStatus;
                sumValues(2) = sumValues(2) + agents2(iAlive).LesionStatus;
                sumValues(3) = sumValues(3) + agents2(iAlive).CancerStatus;
            end
            sumValues(4)  = TestAnnoNEG + TestAnnoHIV;
            sumValues(5)  = median([agents2.Age]);
            sumValues(6)  = safe_div(sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.TestedStatus > 0), ...
                                     sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49));
            sumValues(7)  = safe_div(atnormHIVpos, tottargetHIVpos);
            sumValues(8)  = safe_div(atnormHIVneg, tottargetHIVneg);
            sumValues(9)  = safe_div(sum(agentstable.Gender == 1 & agentstable.Age >= 35 & agentstable.TestedStatus > 1), ...
                                     sum(agentstable.Gender == 1 & agentstable.Age >= 35));
            sumValues(10) = safe_div(incCancer * 100000, denom);
            sumValues(11) = safe_div(sum(agentstable.Gender == 1 & agentstable.InfectionStatus == 1 & agentstable.Age >= 25 & agentstable.Age <= 49), ...
                                     sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49));
            sumValues(12) = safe_div(sum(agentstable.Gender == 1 & agentstable.LesionStatus == 1 & agentstable.Age >= 25 & agentstable.Age <= 49), ...
                                     sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49));
            sumValues(13) = safe_div(incCancerHIV, incCancerNOHIV) * ...
                            safe_div(sum(agentstable.HIVStatus == 0 & agentstable.Gender == 1 & agentstable.Age > 25), ...
                                     sum(agentstable.HIVStatus == 1 & agentstable.Gender == 1 & agentstable.Age > 25));
            sumValues(14) = max(0, safe_div(needed - (TestAnnoNEG + TestAnnoHIV), needed));
            sumValues(15) = safe_div(atnormHIVpos, tottargetHIVpos);
            sumValues(16) = safe_div(atnormHIVneg, tottargetHIVneg);
            sumValues(17) = safe_div(sum(agentstable.Age >= 15 & agentstable.Age <= 49 & agentstable.Gender == 1 & agentstable.HIVStatus > 0), ...
                                     sum(agentstable.Age >= 15 & agentstable.Age <= 49 & agentstable.Gender == 1));
            sumValues(18) = safe_div(sum(agentstable.Gender == 1 & agentstable.Age >= 14 & agentstable.VaccineStatus == 1), ...
                                     sum(agentstable.Gender == 1 & agentstable.Age >= 14));
            sumValues(19) = safe_div(atnormHIVneg + atnormHIVpos, tottargetHIVpos + tottargetHIVneg);
            sumValues(20) = denom;
            sumValues(21) = safe_div(incCancer30 * 100000, sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age < 35));
            sumValues(22) = safe_div(incCancer40 * 100000, sum(agentstable.Gender == 1 & agentstable.Age >= 35 & agentstable.Age < 45));
            sumValues(23) = safe_div(incCancer50 * 100000, sum(agentstable.Gender == 1 & agentstable.Age >= 45 & agentstable.Age < 55));
            sumValues(24) = safe_div(incCancer60 * 100000, sum(agentstable.Gender == 1 & agentstable.Age >= 55));
            lesion0 = safe_div(sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49 & agentstable0.LesionStatus == 1), ...
                               sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49));
            inf0 = safe_div(sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49 & agentstable0.InfectionStatus == 1), ...
                            sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49));
            mat(1,:) = [sum(agentstable0.InfectionStatus), sum(agentstable0.LesionStatus), sum(agentstable0.CancerStatus), 0, ...
                        median(agentstable0.Age), 0, 0, 0, 0, 34, inf0, lesion0, NaN, 1, 1, 1, NaN, 0, 0, populationSize/2, NaN, NaN, NaN, NaN];
            mat(t,:) = sumValues;
            elig(t,1) = length(eligibleagentsHIVpos);
            elig(t,2) = yearlytesttargetHIV(t);
            elig(t,3) = length(eligibleagentsHIVgen);
            elig(t,4) = yearlytesttargetNEG(t);
            elig(t,5) = TestAnnoNEG + TestAnnoHIV;
            I1 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age < 15 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age <= 14));
            I2 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age > 14 & agentstable.Age <= 24 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age > 14 & agentstable.Age <= 24));
            I3 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age > 24 & agentstable.Age <= 34 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age > 24 & agentstable.Age <= 34));
            I4 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age > 34 & agentstable.Age <= 44 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age > 34 & agentstable.Age <= 44));
            I5 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age > 44 & agentstable.Age <= 54 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age > 44 & agentstable.Age <= 54));
            I6 = safe_div(sum(agentstable.Gender == 1 & agentstable.Age > 54 & agentstable.InfectionStatus == 1), ...
                          sum(agentstable.Gender == 1 & agentstable.Age >= 54));
            hpv_by_year(t,:) = [I1, I2, I3, I4, I5, I6];
            popFemaleAge84 = zeros(17,1);
            for aa = 1:17
                lo = 5 * (aa - 1);
                hi = lo + 4;
                popFemaleAge84(aa) = sum(agentstable.Gender == 1 & agentstable.IsAlive == 1 & ...
                                         agentstable.Age >= lo & agentstable.Age <= hi);
            end
            rateAge84 = 100000 * incCancerAge84 ./ max(popFemaleAge84, 1);
            ASR84_by_year(t) = sum(wASR84(:) .* rateAge84(:)) / sum(wASR84);
        end
        % Save all incident cancers from this ensemble.
        cases{nsim} = cases_nsim;
        disp(['Time Step: ' num2str(t)]);
        disp(['Agents with any recorded test: ' num2str(sum([agents.TestedStatus] > 0))]);
        HPVposiINHIVneg = safe_div(sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.InfectionStatus == 1 & agentstable.HIVStatus == 0), ...
                                   sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.HIVStatus == 0));
        HPVposiINHIVpos = safe_div(sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.InfectionStatus == 1 & agentstable.HIVStatus == 1), ...
                                   sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.HIVStatus == 1));
        Tcured(1,:) = [];
        Tcured_real = height(Tcured) - sum(Tcured.CancerStatus);
        disp(['Number of treated individuals: ' num2str(Tcured_real)]);
        M1 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age < 15 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age <= 14));
        M2 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age > 14 & agentstable.Age <= 24 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age > 14 & agentstable.Age <= 24));
        M3 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age > 24 & agentstable.Age <= 34 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age > 24 & agentstable.Age <= 34));
        M4 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age > 34 & agentstable.Age <= 44 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age > 34 & agentstable.Age <= 44));
        M5 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age > 44 & agentstable.Age <= 54 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age > 44 & agentstable.Age <= 54));
        M6 = safe_div(sum(agentstable.Gender == 0 & agentstable.Age > 54 & agentstable.InfectionStatus == 1), ...
                      sum(agentstable.Gender == 0 & agentstable.Age >= 54));
        MPR_HPV = [M1, M2, M3, M4, M5, M6];
        if doPlot
            figure(2)
            subplot(1,2,1)
            plot(hpv_by_year(end,:), 'r')
            hold on
            plot(prevalenceValues0F, 'k')
            legend('Final.Prevalence', 'Initial.Prevalence')
            xticks(1:6)
            xticklabels({'0','14','24','34','44','54'})
            subplot(1,2,2)
            plot(MPR_HPV, 'b')
            hold on
            plot(prevalenceValues0M, 'k')
            legend('Final.Prevalence', 'Initial.Prevalence')
            xticks(1:6)
            xticklabels({'0','14','24','34','44','54'})
        end
        SummaryMat = array2table(mat);
        SummaryMat.Properties.VariableNames(1:24) = {'HPV+','Lesions','Cancer','Test','Age','EverTested', ...
            'Tested.HIVmean','Tested.HIVnegmean','TestedTwice>35','IncCancer','PrevHPV25_49', ...
            'PrevLesion25_49','HPVinHIVratio','MissingTests','CoveragePos','CoverageNeg', ...
            'HIVprev','Vacc_prev','Coverage','Nwomen','Test30','Test40','Test50','Test60'};
        yrs = max(1, timeSteps-9):timeSteps;
        hpv_last10_all(nsim,:) = mean(hpv_by_year(yrs,:), 1, 'omitnan');
        lesion_last10_all(nsim,1) = mean(SummaryMat.PrevLesion25_49(yrs), 'omitnan');
        cancer_last10_all(nsim,1) = mean(ASR84_by_year(yrs), 'omitnan');
        multiplier = realSizeF / sum(agesmoothF);
        CancerPerPop = SummaryMat.IncCancer .* SummaryMat.Nwomen(1:timeSteps) * multiplier / 100000;
        Ntest        = multiplier * elig(:,5);
        ensCancer(:,nsim)       = CancerPerPop;
        ensTested(:,nsim)       = Ntest;
        ensCancer100000(:,nsim) = SummaryMat.IncCancer;
    end
end
% -------------------------------------------------------------------------
function [age19F, age19M] = load_age_profiles(repoRoot)
    candidates = { ...
        fullfile(repoRoot, 'ageprofile.mat'), ...
        fullfile(repoRoot, 'data', 'ageprofile.mat')};
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            S = load(candidates{i});
            if isfield(S, 'age19F') && isfield(S, 'age19M')
                validateattributes(S.age19F, {'numeric'}, ...
                    {'real','finite','nonnegative','vector','numel',17});
                validateattributes(S.age19M, {'numeric'}, ...
                    {'real','finite','nonnegative','vector','numel',17});
                if sum(S.age19F) <= 0 || sum(S.age19M) <= 0
                    error('ZambiaScenarios:InvalidPopulation', ...
                        'age19F and age19M must each have a positive total.');
                end
                age19F = S.age19F(:)';
                age19M = S.age19M(:)';
                return
            end
            error('ageprofile.mat found, but variables age19F/age19M are missing.');
        end
    end
    error('Could not find ageprofile.mat in repo root or data/ subfolder.');
end
function y = safe_div(num, den)
    if den == 0 || isnan(den)
        y = NaN;
    else
        y = num / den;
    end
end
function agents = vaccinate_agents(cover_goal, gender_filter, agentstable, agents)
    AlreadyV = find(gender_filter & agentstable.Age == 13 & agentstable.VaccineStatus == 1);
    eligible_agents_vac = find(gender_filter & agentstable.Age == 13 & agentstable.VaccineStatus == 0);
    AV = numel(AlreadyV);
    EV = numel(eligible_agents_vac);
    pV = AV / max(AV + EV, 1);
    pnow = max(0, cover_goal - pV);
    vyear = round(pnow * (AV + EV));
    if EV == 0 || vyear == 0
        return
    end
    selected_idV = randsample(eligible_agents_vac, min(vyear, EV));
    [agents(selected_idV).VaccineStatus] = deal(1);
end
function band = ageBand84(age)
    if age < 0 || age > 84 || isnan(age)
        band = NaN;
    else
        band = floor(age / 5) + 1;
    end
end
function [agents, agentstable0] = initialization(populationSize, age_vector, agenti, agents)
    for i = 1:populationSize
        currentAge = age_vector(i);
        apippo = agenti(agenti.Age == currentAge, :);
        xx = randperm(size(apippo, 1), 1);
        agents(i).Age = apippo.Age(xx);
        agents(i).Gender = apippo.Gender(xx);
        agents(i).HIVStatus = apippo.HIVStatus(xx);
        agents(i).TestedStatus = apippo.TestedStatus(xx);
        agents(i).LesionStatus = apippo.LesionStatus(xx);
        agents(i).InfectionStatus = apippo.InfectionStatus(xx);
        agents(i).CancerStatus = apippo.CancerStatus(xx);
        agents(i).VaccineStatus = apippo.VaccineStatus(xx);
        agents(i).YearsSinceTested = apippo.YearsSinceTested(xx);
        agents(i).YearsSinceLesions = apippo.YearsSinceLesions(xx);
        agents(i).YearsSinceHPV = apippo.YearsSinceHPV(xx);
        agents(i).IsAlive = apippo.IsAlive(xx);
        agents(i).YearsSinceCancer = apippo.YearsSinceCancer(xx);
        agents(i).TreatmentStatus = 0;
        % Assigned once and retained for the individual's lifetime.
        agents(i).ScreenAccessScore = rand();
    end
    agentstable = struct2table(agents);
    agentstable0 = agentstable(1:populationSize,:);
end
function [population, diff] = growth_rate_paper(initial_population, initial_growth_rate, decrease_rate, time_period)
    population = zeros(time_period, 1);
    diff = zeros(time_period, 1);
    population(1) = initial_population;
    diff(1) = 0;
    grow_rate = initial_growth_rate;
    for t = 2:time_period
        population(t) = population(t-1) * (1 + grow_rate);
        diff(t) = population(t) - population(t-1);
        grow_rate = grow_rate - decrease_rate * grow_rate;
    end
end
function deathProbability2 = storedata()
    WHO_agebandstart = [0 1 5:5:85];
    nMx = [0.04076446
           0.00474964
           0.00130621
           0.00069316
           0.00125168
           0.00191106
           0.00264475
           0.00380903
           0.00521751
           0.00706454
           0.00882847
           0.01154115
           0.01487379
           0.02064697
           0.02920273
           0.04519106
           0.06934874
           0.10908678
           0.19249749];
    deathProbability2 = zeros(1,85);
    for a = 1:85
        ind = max(find([0 1 5:5:85] <= a, 1, 'last'));
        deathProbability2(a) = nMx(ind);
    end
end
