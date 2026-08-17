function [WT_States, w_transition_matrix, siteCaps] = ...
    TwoSite_Copula_Wind_Model( ...
    turbine, WPP_CAP, num_wind_states, ...
    copulaModel, Nsim)

    %% Turbine parameters
    rho_air = 1.225;

    if strcmp(turbine,'IEA')

        wT = 15;

        cutin  = 3;
        cutout = 25;

        Cp = 0.455811125095624;
        Area_wT = pi*120^2;

    elseif strcmp(turbine,'DTU')

        wT = 10;

        cutin  = 4;
        cutout = 25;

        Cp = 0.441353713697127;
        Area_wT = pi*89.15^2;

    else
        error('Invalid turbine type.');
    end

    %% -------------------------------------------------------------
    % Split total installed capacity between the two sites
    % --------------------------------------------------------------

    nT = round(WPP_CAP/wT);

    nT1 = floor(nT/2);
    nT2 = nT - nT1;

    CAP1 = nT1*wT;
    CAP2 = nT2*wT;

    siteCaps = [CAP1 CAP2];

    %% -------------------------------------------------------------
    % Turbine wind-power discrete states
    % --------------------------------------------------------------

    pLevels = wT*linspace(0,1,num_wind_states);       % MW

    wLevels = ((pLevels*1e6)*2 / ...
              (rho_air*Cp*Area_wT)).^(1/3);

    %% -------------------------------------------------------------
    % Generate correlated wind-speed samples from copula
    % -------------------------------------------------------------

    rng(1);          % reproducibility

    U_sim = copularnd( ...
        'Gaussian', ...
        copulaModel.Rho, ...
        Nsim);

    % Recover empirical wind-speed marginal distributions
    w1_sim = empiricalInverse( ...
        copulaModel.w1, U_sim(:,1));

    w2_sim = empiricalInverse( ...
        copulaModel.w2, U_sim(:,2));

    %% -------------------------------------------------------------
    % Convert simulated wind speeds into turbine power states
    % -------------------------------------------------------------

    p1_sim = windToPowerState( ...
        w1_sim, pLevels, wLevels, cutin, cutout);

    p2_sim = windToPowerState( ...
        w2_sim, pLevels, wLevels, cutin, cutout);

    %% Total output assuming all turbines available

    P_sim = nT1*p1_sim + nT2*p2_sim;

    %% -------------------------------------------------------------
    % All possible aggregate wind-output states
    % -------------------------------------------------------------

    pairCaps = nT1*pLevels(:) + nT2*pLevels(:)';

    windCaps = unique(pairCaps(:));

    nStates = length(windCaps);

    %% Copula-derived probabilities

    [~,idx_sim] = ismember(P_sim,windCaps);

    counts = accumarray( ...
        idx_sim,1,[nStates 1]);

    probabilities = counts/sum(counts);

    %% -------------------------------------------------------------
    % Transition rates:
    %
    % Use the ORIGINAL paired time series because copularnd generates
    % samples from the joint distribution but does not preserve the
    % chronological structure of the wind-speed observations.
    % -------------------------------------------------------------

    p1_obs = windToPowerState( ...
        copulaModel.w1, ...
        pLevels,wLevels,cutin,cutout);

    p2_obs = windToPowerState( ...
        copulaModel.w2, ...
        pLevels,wLevels,cutin,cutout);

    P_obs = nT1*p1_obs + nT2*p2_obs;

    [~,idx_obs] = ismember(P_obs,windCaps);

    occurrences = accumarray( ...
        idx_obs,1,[nStates 1]);

    transition_matrix = zeros(nStates);

    from = idx_obs(1:end-1);
    to   = idx_obs(2:end);

    for i = 1:nStates

        if occurrences(i) > 0

            mask = (from == i);

            transitionCounts = accumarray( ...
                to(mask),1,[nStates 1]);

            % Remove self-transitions
            transitionCounts(i) = 0;

            transition_matrix(i,:) = ...
                transitionCounts' / occurrences(i);

        end
    end

    % Assuming the original observations are hourly
    transition_matrix = transition_matrix*24;

    %% -------------------------------------------------------------
    % Remove states never observed/generated
    % -------------------------------------------------------------

    keep = probabilities > 0 | occurrences > 0;

    windCaps = windCaps(keep);
    probabilities = probabilities(keep);

    transition_matrix = ...
        transition_matrix(keep,keep);

    probabilities = probabilities/sum(probabilities);

    %% -------------------------------------------------------------
    % Your original model uses descending capacity ordering
    % -------------------------------------------------------------

    windCaps = flipud(windCaps);
    probabilities = flipud(probabilities);

    w_transition_matrix = ...
        flipud(fliplr(transition_matrix));

    nStates = length(windCaps);

    %% Calculate lambda+ and lambda-

    lambda_plus  = zeros(nStates,1);
    lambda_minus = zeros(nStates,1);

    for i = 1:nStates

        % Higher-capacity states
        if i > 1
            lambda_plus(i) = ...
                sum(w_transition_matrix(i,1:i-1));
        end

        % Lower-capacity states
        if i < nStates
            lambda_minus(i) = ...
                sum(w_transition_matrix(i,i+1:end));
        end

    end

    %% Final aggregate two-WF wind model

    frequency = probabilities .* ...
        (lambda_plus + lambda_minus);

    WT_States = [ ...
        windCaps, ...
        probabilities, ...
        lambda_plus, ...
        lambda_minus, ...
        frequency];

end


%% ========================================================================
function x = empiricalInverse(data,u)

    data = sort(data(:));

    N = length(data);

    q = ((1:N)' - 0.5)/N;

    % Prevent extrapolation beyond empirical distribution
    u = max(u,q(1));
    u = min(u,q(end));

    x = interp1(q,data,u,'linear');

end


%% ========================================================================
function p = windToPowerState( ...
    w,pLevels,wLevels,cutin,cutout)

    p = zeros(size(w));

    active = ...
        w >= cutin & ...
        w <= cutout;

    positiveW = wLevels(2:end);

    if length(positiveW) == 1

        idx = ones(sum(active),1);

    else

        boundaries = ...
            (positiveW(1:end-1) + ...
             positiveW(2:end))/2;

        edges = [-Inf boundaries Inf];

        idx = discretize( ...
            w(active),edges);

    end

    p(active) = pLevels(idx+1);

end