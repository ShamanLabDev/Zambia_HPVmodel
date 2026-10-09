function [ensCancer, ensTested, SummaryMat, ensCancer100000, agentstable, calib] = ...
    ABmodelZambia_paper_calibration_only_repo3(pars, time, ens)
% ABMODELZAMBIA_PAPER_CALIBRATION_ONLY_REPO3 Zambia calibration-only IBM.
%   See README.md for usage, output definitions, and preserved limitations.
%   TIME is 2:81; use 81 for the full 50-row burn-in + 1990:2020 history.
%   The supplied calibrated parameter values and population workbook are
%   required; this function does not estimate parameters.
%   Per-ensemble Twister seed: round(1e6 + 1e3*ensembleIndex + TIME).
%   The caller's RNG state is restored on return or error.
% Required external data file:
%   - Pop_1990.xlsx  (must contain columns named Women and Men)
%
%
% Inputs
%   pars : struct with fields mu, Pcum, a2, a3, a4
%   time : number of yearly simulation steps
%   ens  : number of ensemble simulations
%
% Outputs
%   ensCancer         : yearly incident cancers mapped to real population
%   ensTested         : legacy scaled screening allocation (see README)
%   SummaryMat        : summary table from the final ensemble run
%   ensCancer100000   : yearly cancer incidence per 100,000 women
%   agentstable       : final agent table from the last ensemble run
%   calib             : calibration summaries across ensembles
%
    narginchk(3, 3);
    validateattributes(time, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2, '<=', 81}, mfilename, 'time', 2);
    validateattributes(ens,  {'numeric'}, {'real', 'finite', 'scalar', 'integer', '>=', 1}, mfilename, 'ens', 3);
    validateattributes(pars, {'struct'}, {'scalar'}, mfilename, 'pars', 1);
    requiredFields = {'mu', 'Pcum', 'a2', 'a3', 'a4'};
    for k = 1:numel(requiredFields)
        if ~isfield(pars, requiredFields{k})
            error('ZambiaCalibration:MissingParameter', ...
                'Missing field pars.%s.', requiredFields{k});
        end
    end
    for k = 1:numel(requiredFields)
        validateattributes(pars.(requiredFields{k}), {'numeric'}, ...
            {'real', 'finite', 'scalar', 'nonnegative'}, ...
            mfilename, ['pars.' requiredFields{k}]);
    end
    if pars.mu > 1 || pars.Pcum > 1 || ...
            any([pars.a2, pars.a2 * pars.a3, pars.a2 * pars.a3 * pars.a4] > 1)
        error('ZambiaCalibration:InvalidProbability', ...
            'mu, Pcum, a2, a2*a3, and a2*a3*a4 must be in [0, 1].');
    end

    % Preserve the caller's random stream while retaining the original seeds.
    callerRng = rng;
    restoreRng = onCleanup(@() rng(callerRng));

    burntime = 50;              % first 50 years: fixed age structure
    nens = ens;
    ageGroups = [14, 24, 34, 44, 54];
    timeSteps = time;
    populationSize = 50000;
    % Load 1990 age-sex population profile.
    T = readtable(resolve_data_file('Pop_1990.xlsx'));
    if ~all(ismember({'Women', 'Men'}, T.Properties.VariableNames))
        error('Pop_1990.xlsx must contain columns named Women and Men.');
    end
    age90F = T.Women;
    age90M = T.Men;
    validateattributes(age90F, {'numeric'}, ...
        {'real', 'finite', 'nonnegative', 'column', 'numel', 17}, mfilename, 'Women');
    validateattributes(age90M, {'numeric'}, ...
        {'real', 'finite', 'nonnegative', 'column', 'numel', 17}, mfilename, 'Men');
    if sum(age90F) <= 0 || sum(age90M) <= 0
        error('ZambiaCalibration:InvalidPopulation', ...
            'Women and Men must each have a positive total population.');
    end
    realSizeF = sum(age90F);
    % Historical mortality and total population.
    [mHIV9020, deathProbability2, popZMB_1990_2020] = datasets();
    popZMB_1990_2020(1) = sum(age90F + age90M);
    % Target population by year.
    pop90 = sum(age90F + age90M);
    scalePop = populationSize / pop90;
    normaliz = scalePop;
    populationSizeY = [repmat(populationSize, burntime, 1); round(scalePop * popZMB_1990_2020)];
    maxPopulationSize = max(populationSizeY);
    age90Fyear = smoothdata(repelem(age90F/5, 5), 'movmean', 5);
    age90Myear = smoothdata(repelem(age90M/5, 5), 'movmean', 5);
    agesmoothF = round(age90Fyear(:) * normaliz);
    agesmoothM = round(age90Myear(:) * normaliz);
    % Calibrated parameters.
    probabilityLesion = pars.mu;
    Pcum = pars.Pcum;
    a2 = pars.a2;
    a3 = pars.a3;
    a4 = pars.a4;
    s2 = a2;
    s3 = a2 * a3;
    s4 = a2 * a3 * a4;
    SHPV = [0, s2, s3, s4, 0, 0];
    % Historical HIV incidence anchors (annual probability).
    yrs_anchor  = [1990 1995 2000 2005 2010 2015 2020];
    incF_anchor = [0.039 0.023 0.021 0.0195 0.020 0.010 0.0063];
    incM_anchor = [0.026 0.017 0.0145 0.0125 0.0125 0.0030 0.0006];
    yrs_full  = 1990:2020;
    incF_full = interp1(yrs_anchor, incF_anchor, yrs_full, 'linear');
    incM_full = interp1(yrs_anchor, incM_anchor, yrs_full, 'linear');
    nHistYears = numel(yrs_full);
    % HIV-related multipliers over time.
    RR_acq_start = 2.6;
    RR_acq_end   = 1.7;
    RR_les_start = 3.7;
    RR_les_end   = 2.4;
    HR_clear_start = 0.56;
    HR_clear_end   = 0.70;
    deathProbability = deathProbability2;
    BW = [0, 1.5, 1, 0.8, 0.5, 0.25];
    BM = [0, 1.7, 2, 2, 1.3, 1];
    PartnerFunction = @(Age, Gender) ...
        (Gender == 0) .* ((Age <= ageGroups(1)) * BM(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * BM(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * BM(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * BM(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * BM(5) + ...
        (Age > ageGroups(5)) * BM(6)) + ...
        (Gender == 1) .* ((Age <= ageGroups(1)) * BW(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * BW(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * BW(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * BW(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * BW(5) + ...
        (Age > ageGroups(5)) * BW(6));
    susceptibility = @(Age) (Age <= ageGroups(1)) * SHPV(1) + ...
        (Age > ageGroups(1) & Age <= ageGroups(2)) * SHPV(2) + ...
        (Age > ageGroups(2) & Age <= ageGroups(3)) * SHPV(3) + ...
        (Age > ageGroups(3) & Age <= ageGroups(4)) * SHPV(4) + ...
        (Age > ageGroups(4) & Age <= ageGroups(5)) * SHPV(5) + ...
        (Age > ageGroups(5)) * SHPV(6);
    recoveryrate1 = @(YearsSinceHPV) (YearsSinceHPV == 1) * 0.6 + ...
                                    (YearsSinceHPV == 2) * 0.5 + ...
                                    (YearsSinceHPV == 3) * 0.2 + ...
                                    (YearsSinceHPV == 4) * 0.2;
    Py = 1 - (1 - Pcum)^(1/15);
    CancerRatefromLesion = @(YearsSinceLesions) (YearsSinceLesions > 5) * Py;
    DeathCancer = @(YearsSinceCancer) (YearsSinceCancer <= 5) * 0.12 + ...
                                      (YearsSinceCancer <= 10 & YearsSinceCancer > 5) * 0.05;
    wASR84 = [12000 10000 9000 9000 8000 8000 6000 6000 6000 ...
              6000 5000 4000 4000 3000 2000 1000 500];
    multi = sum(age90F, 1) / sum(agesmoothF);
    hpv_last5_all    = zeros(nens, 6);
    lesion_last5_all = zeros(nens, 1);
    cancer_last5_all = zeros(nens, 1);
    % Historical screening allocation (VIA, 2006:2020).
    sensitivity = 0.35;
    ratiotreated = 0.73;
    efficacyHIV = 0.63;
    efficacyNEG = 0.85;
    probabilityCuredHIV = sensitivity * efficacyHIV * ratiotreated;
    probabilityCuredNEG = sensitivity * efficacyNEG * ratiotreated;
    tin = burntime + 16;
    screen_total = [zeros(1, tin), 8000 9500 11000 13000 15352 16000 15000 15090 18000 20000 22000 24000 26000 32783 32783];
    screen_hiv   = [zeros(1, tin), 2480 2945 3410 4030 5757 6000 5625 5659 6750 7500 8250 9000 9750 12294 12294];
    screen_else  = screen_total - screen_hiv;
    ensCancer = zeros(timeSteps, nens);
    ensTested = zeros(timeSteps, nens);
    ensCancer100000 = zeros(timeSteps, nens);
    % Retain annual calibration diagnostics for every ensemble simulation.
    cancerASRByYear_all = NaN(timeSteps, nens);
    hpvByYear_all       = NaN(timeSteps, 6, nens);
    lesionByYear_all    = NaN(timeSteps, nens);
    for nsim = 1:nens
        baseSeed = round(1e6 + 1e3 * nsim + timeSteps);
        rng(baseSeed, 'twister');
        % Reset annual diagnostic trajectory for this ensemble simulation.
        ASR84_by_year = NaN(timeSteps, 1);
        emptyAgent = struct(...
            'Age', NaN, ...
            'Gender', NaN, ...
            'HIVStatus', NaN, ...
            'VaccineStatus', NaN, ...
            'InfectionStatus', NaN, ...
            'TestedStatus', NaN, ...
            'LesionStatus', NaN, ...
            'CancerStatus', NaN, ...
            'YearsSinceTested', NaN, ...
            'IsAlive', NaN, ...
            'YearsSinceHPV', NaN, ...
            'YearsSinceLesions', NaN, ...
            'YearsSinceCancer', NaN, ...
            'TreatmentStatus', NaN);
        agents = repmat(emptyAgent, maxPopulationSize, 1);
        [agents, agentstable0] = initializeagents(agesmoothF, agesmoothM, populationSize, agents);
        targetAgeSex = zeros(85, 2);
        for a = 0:84
            for g = 0:1
                targetAgeSex(a+1, g+1) = sum(agentstable0.IsAlive == 1 & agentstable0.Age == a & agentstable0.Gender == g);
            end
        end
        baselinePool = cell(85, 2);
        for a = 0:84
            for g = 0:1
                baselinePool{a+1, g+1} = agentstable0(agentstable0.IsAlive == 1 & agentstable0.Age == a & agentstable0.Gender == g, :);
            end
        end
        mat = zeros(timeSteps, 24);
        hpv_by_year = nan(timeSteps, 6);
        Ncured = 0;
        female0 = agentstable0(agentstable0.Gender == 1 & agentstable0.IsAlive == 1, :);
        hpv_by_year(1,1) = mean(female0.InfectionStatus(female0.Age <= 14));
        hpv_by_year(1,2) = mean(female0.InfectionStatus(female0.Age >= 15 & female0.Age <= 24));
        hpv_by_year(1,3) = mean(female0.InfectionStatus(female0.Age >= 25 & female0.Age <= 34));
        hpv_by_year(1,4) = mean(female0.InfectionStatus(female0.Age >= 35 & female0.Age <= 44));
        hpv_by_year(1,5) = mean(female0.InfectionStatus(female0.Age >= 45 & female0.Age <= 54));
        hpv_by_year(1,6) = mean(female0.InfectionStatus(female0.Age >= 55));
        hpv_by_year(1, isnan(hpv_by_year(1,:))) = 0;
        Tcancer = table();
        Tcured = table();
        cancerDeathCounter = 0;
        elig = zeros(timeSteps, 5);
        elig(1,:) = 0;
        for t = 2:timeSteps
            incCancerAge84 = zeros(17, 1);
            incCancer = 0;
            incCancerHIV = 0;
            incCancerNOHIV = 0;
            incCancer20 = 0;
            incCancer30 = 0;
            incCancer40 = 0;
            incCancer50 = 0;
            incCancer60 = 0;
            incCancerDeath = 0;
            agentstable = struct2table(agents);
            state0 = agents;
            nsexvec = poissrnd(PartnerFunction(agentstable.Age, agentstable.Gender));
            if t > burntime
                w = (t - burntime) / (timeSteps - burntime);
                histIdx = min(t - burntime, nHistYears);
                incF_t = incF_full(histIdx);
                incM_t = incM_full(histIdx);
                mHIV = mHIV9020(histIdx);
            else
                w = 0;
                incF_t = mHIV9020(1) * 0.3;
                incM_t = mHIV9020(1) * 0.3;
                mHIV = mHIV9020(1);
            end
            incHIVM = [0, incM_t, incM_t, incM_t, 0, 0];
            incHIVF = [0, incF_t, incF_t, incF_t, 0, 0];
            HIVincFunction = @(Age, Gender) ...
                (Gender == 0) .* ((Age <= ageGroups(1)) * incHIVM(1) + ...
                (Age > ageGroups(1) & Age <= ageGroups(2)) * incHIVM(2) + ...
                (Age > ageGroups(2) & Age <= ageGroups(3)) * incHIVM(3) + ...
                (Age > ageGroups(3) & Age <= ageGroups(4)) * incHIVM(4) + ...
                (Age > ageGroups(4) & Age <= ageGroups(5)) * incHIVM(5) + ...
                (Age > ageGroups(5)) * incHIVM(6)) + ...
                (Gender == 1) .* ((Age <= ageGroups(1)) * incHIVF(1) + ...
                (Age > ageGroups(1) & Age <= ageGroups(2)) * incHIVF(2) + ...
                (Age > ageGroups(2) & Age <= ageGroups(3)) * incHIVF(3) + ...
                (Age > ageGroups(3) & Age <= ageGroups(4)) * incHIVF(4) + ...
                (Age > ageGroups(4) & Age <= ageGroups(5)) * incHIVF(5) + ...
                (Age > ageGroups(5)) * incHIVF(6));
            RR_acq = exp((1 - w) * log(RR_acq_start) + w * log(RR_acq_end));
            RR_lesion = exp((1 - w) * log(RR_les_start) + w * log(RR_les_end));
            HR_clear = exp((1 - w) * log(HR_clear_start) + w * log(HR_clear_end));
            for i = 1:maxPopulationSize
                if isnan(agents(i).IsAlive) || agents(i).IsAlive ~= 1
                    continue
                end
                agents(i).Age = agents(i).Age + 1;
                if agents(i).Age >= 85 || agents(i).IsAlive == 0
                    agents(i).IsAlive = 0;
                    continue
                end
                age = agents(i).Age;
                deathP = deathProbability(max(round(age), 1)) + mHIV * agents(i).HIVStatus;
                if rand() < deathP
                    agents(i).IsAlive = 0;
                    continue
                end
                if agents(i).Age >= 15
                    if rand() < HIVincFunction(agents(i).Age, agents(i).Gender) + agents(i).HIVStatus
                        agents(i).HIVStatus = 1;
                    end
                    if agents(i).Gender == 0
                        if agents(i).InfectionStatus == 0
                            nsex = nsexvec(i);
                            if nsex > 0
                                eligible_agents_idx = find(agentstable.Gender == 1 & (age - agentstable.Age) <= 5 & (age - agentstable.Age) >= -1 & agentstable.Age >= 15);
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
                                        agents(i).YearsSinceHPV = 1;
                                    end
                                end
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 0
                            transitionDraw = rand();
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
                            if transitionDraw < p1
                                agents(i).LesionStatus = 1;
                                agents(i).YearsSinceLesions = 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            elseif transitionDraw < p1 + pclr
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1
                            transitionDraw = rand();
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if transitionDraw < pclr
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
                        if agents(i).InfectionStatus == 0
                            nsex = nsexvec(i);
                            if nsex > 0
                                eligible_agents_idx = find(agentstable.Gender == 0 & (agentstable.Age - age) <= 5 & (agentstable.Age - age) >= -1 & agentstable.Age >= 15);
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
                                        agents(i).YearsSinceHPV = 1;
                                    end
                                end
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 0
                            if agents(i).YearsSinceTested == 0 && agents(i).YearsSinceHPV >= 1
                                agents(i).YearsSinceTested = 6;
                            end
                            transitionDraw = rand();
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
                            if transitionDraw < p1
                                agents(i).LesionStatus = 1;
                                agents(i).YearsSinceLesions = 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            elseif transitionDraw < p1 + pclr
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1 && agents(i).CancerStatus == 0
                            transitionDraw = rand();
                            p2 = CancerRatefromLesion(agents(i).YearsSinceLesions);
                            pclr0 = recoveryrate1(agents(i).YearsSinceHPV);
                            if agents(i).HIVStatus == 1
                                pclr = 1 - (1 - pclr0)^HR_clear;
                            else
                                pclr = pclr0;
                            end
                            if transitionDraw < p2
                                agents(i).CancerStatus = 1;
                                incCancer = incCancer + 1;
                                incCancerHIV = incCancerHIV + agents(i).HIVStatus;
                                incCancerNOHIV = incCancerNOHIV + (1 - agents(i).HIVStatus);
                                band = ageBand84(agents(i).Age);
                                if ~isnan(band)
                                    incCancerAge84(band) = incCancerAge84(band) + 1;
                                end
                                incCancer20 = incCancer20 + (agents(i).Age >= 15 & agents(i).Age < 25);
                                incCancer30 = incCancer30 + (agents(i).Age >= 25 & agents(i).Age < 35);
                                incCancer40 = incCancer40 + (agents(i).Age >= 35 & agents(i).Age < 45);
                                incCancer50 = incCancer50 + (agents(i).Age >= 45 & agents(i).Age < 55);
                                incCancer60 = incCancer60 + (agents(i).Age >= 55);
                                agents(i).YearsSinceCancer = 1;
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                            elseif transitionDraw < p2 + pclr
                                agents(i).LesionStatus = 0;
                                agents(i).YearsSinceLesions = 0;
                                agents(i).InfectionStatus = 0;
                                agents(i).YearsSinceHPV = 0;
                            else
                                agents(i).YearsSinceLesions = agents(i).YearsSinceLesions + 1;
                                agents(i).YearsSinceHPV = agents(i).YearsSinceHPV + 1;
                            end
                        elseif agents(i).InfectionStatus == 1 && agents(i).LesionStatus == 1 && agents(i).CancerStatus == 1
                            transitionDraw = rand();
                            p3 = DeathCancer(agents(i).YearsSinceCancer) * (1 + agents(i).HIVStatus * 0.1);
                            if transitionDraw < p3
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
            if t <= burntime
                agents = rebalance_age_structure_exact_age(agents, targetAgeSex, baselinePool);
            else
                agents = fill_to_target_with_newborns(agents, populationSizeY(t));
            end
            Positive_after_treatment = 0.1;
            agentstable = struct2table(agents);
            if t > timeSteps - 14 && t <= timeSteps
                eligibleagentsHIVgen = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 5 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 0);
                eligibleagentsHIVpos = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 3 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 1);
                needed = length(eligibleagentsHIVpos) + length(eligibleagentsHIVgen);
                TestAnnoHIV = min(round(screen_hiv(t) * 0.0068), length(eligibleagentsHIVpos));
                TestAnnoNEG = min(round(screen_else(t) * 0.0068), length(eligibleagentsHIVgen));
                testwinnersID_HIV = randperm(length(eligibleagentsHIVpos), TestAnnoHIV);
                testwinnersID_NEG = randperm(length(eligibleagentsHIVgen), TestAnnoNEG);
                for ii = testwinnersID_HIV
                    yy = eligibleagentsHIVpos(ii);
                    agents(yy).TestedStatus = agents(yy).TestedStatus + 1;
                    agents(yy).YearsSinceTested = 0;
                    if agents(yy).LesionStatus == 1 && agents(yy).CancerStatus == 0
                        agents(yy).TreatmentStatus = 1;
                        if rand() < probabilityCuredHIV
                            agents(yy).LesionStatus = 0;
                            agents(yy).YearsSinceLesions = 0;
                            Ncured = Ncured + 1;
                            Tcured = vertcat(Tcured, agentstable(yy,:));
                            if rand() < (1 - Positive_after_treatment)
                                agents(yy).InfectionStatus = 0;
                            end
                        end
                    end
                end
                for ii = testwinnersID_NEG
                    yy = eligibleagentsHIVgen(ii);
                    agents(yy).TestedStatus = agents(yy).TestedStatus + 1;
                    agents(yy).YearsSinceTested = 0;
                    if agents(yy).LesionStatus == 1 && agents(yy).CancerStatus == 0
                        agents(yy).TreatmentStatus = 1;
                        if rand() < probabilityCuredNEG
                            agents(yy).LesionStatus = 0;
                            agents(yy).YearsSinceLesions = 0;
                            Ncured = Ncured + 1;
                            Tcured = vertcat(Tcured, agentstable(yy,:));
                            if rand() < (1 - Positive_after_treatment)
                                agents(yy).InfectionStatus = 0;
                            end
                        end
                    end
                end
                agentstable = struct2table(agents);
                eligibleagentsHIVgen = find(agentstable.Gender == 1 & agentstable.YearsSinceTested > 5 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1);
            else
                TestAnnoNEG = 0;
                TestAnnoHIV = 0;
                eligibleagentsHIVpos = 0;
                eligibleagentsHIVgen = 0;
                needed = 0;
            end
            atnormHIVgen = sum(agentstable.Gender == 1 & agentstable.YearsSinceTested <= 5 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 0);
            atnormHIVpos = sum(agentstable.Gender == 1 & agentstable.YearsSinceTested <= 3 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 1);
            tottargetHIVneg = sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 0);
            tottargetHIVpos = sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.IsAlive == 1 & agentstable.HIVStatus == 1);
            sumValues = zeros(1, 24);
            agents2 = agents([agents.IsAlive] == 1 & [agents.Gender] == 1);
            denom = sum(agentstable.Gender == 1);
            for i = 1:length(agents2)
                sumValues(1) = sumValues(1) + agents2(i).InfectionStatus;
                sumValues(2) = sumValues(2) + agents2(i).LesionStatus;
                sumValues(3) = sumValues(3) + agents2(i).CancerStatus;
            end
            sumValues(4) = TestAnnoNEG + TestAnnoHIV;
            sumValues(5) = median([agents2.Age]);
            sumValues(6) = sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49 & agentstable.TestedStatus > 0) / sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49);
            sumValues(7) = atnormHIVpos / max(tottargetHIVpos, 1);
            sumValues(8) = atnormHIVgen / max(tottargetHIVneg, 1);
            sumValues(9) = sum(agentstable.Gender == 1 & agentstable.Age >= 35 & agentstable.TestedStatus > 1) / sum(agentstable.Gender == 1 & agentstable.Age >= 35);
            sumValues(10) = incCancer * 100000 / max(denom, 1);
            sumValues(11) = sum(agentstable.Gender == 1 & agentstable.InfectionStatus == 1 & agentstable.Age >= 25 & agentstable.Age <= 49) / sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49);
            sumValues(12) = sum(agentstable.Gender == 1 & agentstable.LesionStatus == 1 & agentstable.Age >= 25 & agentstable.Age <= 49) / sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age <= 49);
            if incCancerNOHIV > 0 && sum(agentstable.HIVStatus == 1 & agentstable.Gender == 1 & agentstable.Age > 25) > 0
                sumValues(13) = incCancerHIV / incCancerNOHIV * sum(agentstable.HIVStatus == 0 & agentstable.Gender == 1 & agentstable.Age > 25) / sum(agentstable.HIVStatus == 1 & agentstable.Gender == 1 & agentstable.Age > 25);
            else
                sumValues(13) = NaN;
            end
            if needed > 0
                sumValues(14) = max(0, (needed - (TestAnnoNEG + TestAnnoHIV)) / needed);
            else
                sumValues(14) = 0;
            end
            sumValues(15) = atnormHIVpos / max(tottargetHIVpos, 1);
            sumValues(16) = atnormHIVgen / max(tottargetHIVneg, 1);
            sumValues(17) = sum(agentstable.Age >= 15 & agentstable.Age <= 49 & agentstable.Gender == 1 & agentstable.HIVStatus > 0) / sum(agentstable.Age >= 15 & agentstable.Age <= 49 & agentstable.Gender == 1);
            sumValues(18) = sum(agentstable.Gender == 1 & agentstable.Age >= 14 & agentstable.VaccineStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age >= 14);
            sumValues(19) = (atnormHIVgen + atnormHIVpos) / max(tottargetHIVpos + tottargetHIVneg, 1);
            sumValues(20) = denom;
            sumValues(21) = incCancer30 * 100000 / max(sum(agentstable.Gender == 1 & agentstable.Age >= 25 & agentstable.Age < 35), 1);
            sumValues(22) = incCancer40 * 100000 / max(sum(agentstable.Gender == 1 & agentstable.Age >= 35 & agentstable.Age < 45), 1);
            sumValues(23) = incCancer50 * 100000 / max(sum(agentstable.Gender == 1 & agentstable.Age >= 45 & agentstable.Age < 55), 1);
            sumValues(24) = incCancer60 * 100000 / max(sum(agentstable.Gender == 1 & agentstable.Age >= 55), 1);
            lesion0 = sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49 & agentstable0.LesionStatus == 1) / sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49);
            inf0 = sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49 & agentstable0.InfectionStatus == 1) / sum(agentstable0.Gender == 1 & agentstable0.Age >= 25 & agentstable0.Age <= 49);
            mat(1,:) = [sum(agentstable0.InfectionStatus), sum(agentstable0.LesionStatus), sum(agentstable0.CancerStatus), 0, ...
                        median(agentstable0.Age), 0, 0, 0, 0, 34, inf0, lesion0, NaN, 1, 1, 1, NaN, 0, 0, populationSize/2, NaN, NaN, NaN, NaN];
            mat(t,:) = sumValues;
            elig(t,1) = length(eligibleagentsHIVpos);
            elig(t,2) = screen_hiv(t);
            elig(t,3) = length(eligibleagentsHIVgen);
            elig(t,4) = screen_else(t);
            elig(t,5) = screen_total(t);
            I1 = sum(agentstable.Gender == 1 & agentstable.Age < 15 & agentstable.InfectionStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age <= 14);
            I2 = sum(agentstable.Gender == 1 & agentstable.Age > 14 & agentstable.Age <= 24 & agentstable.InfectionStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age > 14 & agentstable.Age <= 24);
            I3 = sum(agentstable.Gender == 1 & agentstable.Age > 24 & agentstable.Age <= 34 & agentstable.InfectionStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age > 24 & agentstable.Age <= 34);
            I4 = sum(agentstable.Gender == 1 & agentstable.Age > 34 & agentstable.Age <= 44 & agentstable.InfectionStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age > 34 & agentstable.Age <= 44);
            I5 = sum(agentstable.Gender == 1 & agentstable.Age > 44 & agentstable.Age <= 54 & agentstable.InfectionStatus == 1) / sum(agentstable.Gender == 1 & agentstable.Age > 44 & agentstable.Age <= 54);
            I6 = sum(agentstable.Gender == 1 & agentstable.Age >= 55 & agentstable.InfectionStatus == 1) / ...
                 sum(agentstable.Gender == 1 & agentstable.Age >= 55);
            hpv_by_year(t,:) = [I1, I2, I3, I4, I5, I6];
            popFemaleAge84 = zeros(17,1);
            for aa = 1:17
                lo = 5*(aa-1);
                hi = lo + 4;
                popFemaleAge84(aa) = sum(agentstable.Gender == 1 & agentstable.IsAlive == 1 & agentstable.Age >= lo & agentstable.Age <= hi);
            end
            rateAge84 = 100000 * incCancerAge84 ./ max(popFemaleAge84, 1);
            ASR84_by_year(t) = sum(wASR84(:) .* rateAge84(:)) / sum(wASR84);
        end
        cancerASRByYear_all(:, nsim) = ASR84_by_year(:);
        hpvByYear_all(:, :, nsim)    = hpv_by_year;
        disp(['Time Step: ' num2str(t)]);
        disp(['Agents with any recorded test: ' num2str(sum([agents.TestedStatus] > 0))]);
        SummaryMat = array2table(mat);
        SummaryMat.Properties.VariableNames(1:24) = {'HPV+','Lesions','Cancer','Test','Age','EverTested','Tested.HIVmean','Tested.HIVnegmean','TestedTwice>35','IncCancer','PrevHPV25_49', ...
            'PrevLesion25_49','HPVinHIVratio','MissingTests','CoveragePos','CoverageNeg','HIVprev','Vacc_prev','Coverage','Nwomen','T30','T40','T50','T60'};
        % Save annual lesion prevalence for this ensemble simulation.
        lesionByYear_all(:, nsim) = SummaryMat.PrevLesion25_49(:);
        yrs = max(1, timeSteps-4):timeSteps;
        hpv_last5_all(nsim,:)    = mean(hpv_by_year(yrs,:), 1, 'omitnan');
        lesion_last5_all(nsim,1) = mean(SummaryMat.PrevLesion25_49(yrs), 'omitnan');
        cancer_last5_all(nsim,1) = mean(ASR84_by_year(yrs), 'omitnan');
        multiplier = realSizeF / sum(agesmoothF);
        CancerPerPop = SummaryMat.IncCancer .* SummaryMat.Nwomen(1:timeSteps) * multiplier / 100000;
        Ntest = multiplier * elig(:,5);
        ensCancer(:,nsim) = CancerPerPop;
        ensTested(:,nsim) = Ntest;
        ensCancer100000(:,nsim) = SummaryMat.IncCancer;
    end
    calib.hpv    = mean(hpv_last5_all, 1, 'omitnan');
    calib.lesion = mean(lesion_last5_all, 1, 'omitnan');
    calib.cancer = mean(cancer_last5_all, 1, 'omitnan');
    calib.hpv_by_ens    = hpv_last5_all;
    calib.lesion_by_ens = lesion_last5_all;
    calib.cancer_by_ens = cancer_last5_all;
    calib.cancerASRByYear = cancerASRByYear_all;
    calib.hpvByYear       = hpvByYear_all;
    calib.lesionByYear    = lesionByYear_all;
end
function fullpath = resolve_data_file(filename)
    here = fileparts(mfilename('fullpath'));
    candidates = {fullfile(here, filename), fullfile(here, 'data', filename)};
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file')
            fullpath = candidates{i};
            return
        end
    end
    error('ZambiaCalibration:MissingData', ...
        'Required data file %s must be beside this function or in its data folder.', filename);
end
function band = ageBand84(age)
    if age < 0 || age > 84 || isnan(age)
        band = NaN;
    else
        band = floor(age/5) + 1;
    end
end
function agents = rebalance_age_structure_exact_age(agents, targetAgeSex, baselinePool)
    T = struct2table(agents);
    aliveMask = T.IsAlive == 1;
    reassignIdx = find(~aliveMask);
    for a = 0:84
        for g = 0:1
            idx = find(aliveMask & T.Age == a & T.Gender == g);
            nNow = numel(idx);
            nTar = targetAgeSex(a+1, g+1);
            if nNow > nTar
                extra = nNow - nTar;
                killLocal = idx(randperm(nNow, extra));
                aliveMask(killLocal) = false;
                [agents(killLocal).IsAlive] = deal(0);
                reassignIdx = [reassignIdx; killLocal(:)]; %#ok<AGROW>
            end
        end
    end
    T = struct2table(agents);
    T.IsAlive = aliveMask;
    slotPtr = 1;
    for a = 0:84
        for g = 0:1
            idx = find(aliveMask & T.Age == a & T.Gender == g);
            nNow = numel(idx);
            nTar = targetAgeSex(a+1, g+1);
            if nNow < nTar
                need = nTar - nNow;
                for k = 1:need
                    if slotPtr > numel(reassignIdx)
                        error('Not enough slots available in rebalance_age_structure_exact_age');
                    end
                    ii = reassignIdx(slotPtr);
                    slotPtr = slotPtr + 1;
                    donorRow = draw_same_age_donor(T, aliveMask, a, g, baselinePool);
                    agents(ii).Age = a;
                    agents(ii).Gender = g;
                    agents(ii).HIVStatus = donorRow.HIVStatus;
                    agents(ii).VaccineStatus = donorRow.VaccineStatus;
                    agents(ii).InfectionStatus = donorRow.InfectionStatus;
                    agents(ii).TestedStatus = donorRow.TestedStatus;
                    agents(ii).LesionStatus = donorRow.LesionStatus;
                    agents(ii).CancerStatus = donorRow.CancerStatus;
                    agents(ii).YearsSinceTested = donorRow.YearsSinceTested;
                    agents(ii).IsAlive = 1;
                    agents(ii).YearsSinceHPV = donorRow.YearsSinceHPV;
                    agents(ii).YearsSinceLesions = donorRow.YearsSinceLesions;
                    agents(ii).YearsSinceCancer = donorRow.YearsSinceCancer;
                    agents(ii).TreatmentStatus = 0;
                    aliveMask(ii) = true;
                    T.Age(ii) = agents(ii).Age;
                    T.Gender(ii) = agents(ii).Gender;
                    T.HIVStatus(ii) = agents(ii).HIVStatus;
                    T.VaccineStatus(ii) = agents(ii).VaccineStatus;
                    T.InfectionStatus(ii) = agents(ii).InfectionStatus;
                    T.TestedStatus(ii) = agents(ii).TestedStatus;
                    T.LesionStatus(ii) = agents(ii).LesionStatus;
                    T.CancerStatus(ii) = agents(ii).CancerStatus;
                    T.YearsSinceTested(ii) = agents(ii).YearsSinceTested;
                    T.IsAlive(ii) = 1;
                    T.YearsSinceHPV(ii) = agents(ii).YearsSinceHPV;
                    T.YearsSinceLesions(ii) = agents(ii).YearsSinceLesions;
                    T.YearsSinceCancer(ii) = agents(ii).YearsSinceCancer;
                    T.TreatmentStatus(ii) = agents(ii).TreatmentStatus;
                end
            end
        end
    end
    notAlive = ~aliveMask;
    [agents(notAlive).IsAlive] = deal(0);
end
function donorRow = draw_same_age_donor(T, aliveMask, a, g, baselinePool)
    cand = find(aliveMask & T.Age == a & T.Gender == g);
    if ~isempty(cand)
        jj = cand(randi(numel(cand)));
        donorRow = T(jj,:);
        return
    end
    B = baselinePool{a+1, g+1};
    if ~isempty(B)
        jj = randi(height(B));
        donorRow = B(jj,:);
        return
    end
    error('No exact-age donor available for age=%d sex=%d', a, g);
end
function agents = fill_to_target_with_newborns(agents, targetPop)
    T = struct2table(agents);
    aliveIdx = find(T.IsAlive == 1);
    nAlive = numel(aliveIdx);
    if nAlive > targetPop
        warning('Alive population (%d) already exceeds targetPop (%d). No newborns added this year.', nAlive, targetPop);
        return
    end
    need = targetPop - nAlive;
    if need == 0
        return
    end
    freeIdx = find(isnan(T.IsAlive) | T.IsAlive == 0);
    if numel(freeIdx) < need
        error('Not enough free slots to add newborns: need %d, have %d', need, numel(freeIdx));
    end
    fillIdx = freeIdx(1:need);
    for j = 1:numel(fillIdx)
        ii = fillIdx(j);
        agents(ii).Age = 0;
        agents(ii).Gender = double(rand() < 0.5);
        agents(ii).HIVStatus = 0;
        agents(ii).VaccineStatus = 0;
        agents(ii).InfectionStatus = 0;
        agents(ii).TestedStatus = 0;
        agents(ii).LesionStatus = 0;
        agents(ii).CancerStatus = 0;
        agents(ii).YearsSinceTested = 100;
        agents(ii).IsAlive = 1;
        agents(ii).YearsSinceHPV = 0;
        agents(ii).YearsSinceLesions = 0;
        agents(ii).YearsSinceCancer = 0;
        agents(ii).TreatmentStatus = 0;
    end
end
function b = ageBand6(a)
    if a <= 14
        b = 1;
    elseif a <= 24
        b = 2;
    elseif a <= 34
        b = 3;
    elseif a <= 44
        b = 4;
    elseif a <= 54
        b = 5;
    else
        b = 6;
    end
end
function [agents, agentstable0] = initializeagents(agesmoothF, agesmoothM, populationSize, agents)
    age = (0:84)';
    ageF0 = repelem(age, agesmoothF(:));
    ageF = [ageF0; 0; 0; 0];
    ageM = repelem(age, agesmoothM(:));
    AgeAll = [ageF; ageM];
    if numel(AgeAll) < populationSize
        error('ZambiaCalibration:InsufficientInitialPopulation', ...
            ['Smoothed population plus the original three newborns has %d agents; ' ...
             '%d are required. Check the input profile and rounding.'], ...
             numel(AgeAll), populationSize);
    end
    GenderAll = [ones(numel(ageF),1); zeros(numel(ageM),1)];
    perm = randperm(numel(AgeAll));
    AgeAll = AgeAll(perm);
    GenderAll = GenderAll(perm);
    hpv0F = [0.00 0.35 0.30 0.20 0.15 0.15]*1.15;
    hpv0M = [0.00 0.30 0.40 0.35 0.25 0.25]*1.15;
    hiv0F = [0.00 0.03 0.10 0.12 0.09 0.05];
    hiv0M = [0.00 0.01 0.05 0.07 0.05 0.03];
    for i = 1:populationSize
        a = AgeAll(i);
        g = GenderAll(i);
        b = ageBand6(a);
        agents(i).Age = a;
        agents(i).Gender = g;
        agents(i).IsAlive = 1;
        agents(i).VaccineStatus = 0;
        agents(i).TreatmentStatus = 0;
        agents(i).TestedStatus = 0;
        agents(i).YearsSinceTested = 100;
        if g == 1
            pHIV = hiv0F(b);
            pHPV = hpv0F(b);
        else
            pHIV = hiv0M(b);
            pHPV = hpv0M(b);
        end
        agents(i).HIVStatus = rand < pHIV;
        pHPV = min(pHPV * (1 + 0.35 * agents(i).HIVStatus), 0.95);
        if a < 15
            pHPV = 0;
        end
        agents(i).InfectionStatus = rand < pHPV;
        if agents(i).InfectionStatus
            if g == 1
                pLes = 0.02 + 0.05*(a >= 25) + 0.05*agents(i).HIVStatus;
                pLes = min(pLes, 0.35);
            else
                pLes = 0.01 + 0.03*agents(i).HIVStatus;
                pLes = min(pLes, 0.15);
            end
            agents(i).LesionStatus = rand < pLes;
            agents(i).YearsSinceHPV = randi([1 4]);
        else
            agents(i).LesionStatus = 0;
            agents(i).YearsSinceHPV = 0;
        end
        agents(i).CancerStatus = 0;
        if g == 1 && agents(i).LesionStatus && a >= 35
            pCan = 0.001*(a < 45) + 0.003*(a >= 45 && a < 55) + 0.008*(a >= 55);
            agents(i).CancerStatus = rand < pCan;
        end
        if agents(i).LesionStatus
            agents(i).YearsSinceLesions = randi([1 8]);
        else
            agents(i).YearsSinceLesions = 0;
        end
        if agents(i).CancerStatus
            agents(i).InfectionStatus = 1;
            agents(i).LesionStatus = 1;
            agents(i).YearsSinceHPV = max(agents(i).YearsSinceHPV, randi([3 10]));
            agents(i).YearsSinceLesions = max(agents(i).YearsSinceLesions, randi([3 12]));
            agents(i).YearsSinceCancer = randi([1 5]);
        else
            agents(i).YearsSinceCancer = 0;
        end
    end
    agentstable0 = struct2table(agents(1:populationSize));
end
function [mHIV9020, deathProbability2, popZMB_1990_2020] = datasets()
    WHO_agebandstart = [0 1 5:5:85];
    mHIV9020 = [0.024535785
    0.027410113
    0.030665764
    0.034283698
    0.038095475
    0.042046913
    0.046000826
    0.049955598
    0.053914452
    0.057587694
    0.061039698
    0.064107335
    0.06721595
    0.06563013
    0.061051599
    0.056195643
    0.050232475
    0.038631565
    0.026458858
    0.020907143
    0.019374044
    0.017832299
    0.015685955
    0.015245017
    0.015416969
    0.014858073
    0.013330555
    0.012837653
    0.013533899
    0.013427414
    0.013006309
    0.01178026];
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
    nqx = [0.19249749
    0.03963351
    0.01878445
    0.00650980
    0.00345978
    0.00623887
    0.00950984
    0.01313687
    0.01886550
    0.02575166
    0.03470968
    0.04318913
    0.05608749
    0.07170272
    0.09816770
    0.13607896
    0.20301872
    0.29551049
    0.42855869
    1.00000000];
    ex = [65.37
    67.06
    64.31
    59.71
    54.91
    50.24
    45.70
    41.27
    37.02
    32.93
    29.03
    25.22
    21.57
    18.05
    14.74
    11.67
    9.00
    6.73
    4.91];
    deathProbability2 = zeros(1, 85);
    for a = 1:85
        ind = max(find(WHO_agebandstart <= a, 1, 'last'));
        deathProbability2(a) = nMx(ind);
    end
    popZMB_1990_2020 = [7786169
        7981650
        8176680
        8373921
        8576269
        8785763
        9004053
        9237063
        9482408
        9740005
        10017631
        10325185
        10647949
        10983595
        11338198
        11718819
        12129553
        12565085
        13021324
        13490389
        13965594
        14437796
        14913629
        15398997
        15895315
        16399089
        16914423
        17441320
        17973569
        18513839
        19059395];
end
