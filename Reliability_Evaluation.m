function [LOLP, LOLE, LOEE, EIR, Wind_Eshare, EENS_0, EENS_W] = Reliability_Evaluation(load_model, SCOPT_trimmed, Annual_peak, CONVOLVED_SCOPT_merged, SCOPT)
% Reliability_Evaluation calculates various reliability indices for power systems.
% This function supports two load modeling approaches: 'cumulative' and 'exact'.
% It computes Loss of Load Probability (LOLP), Loss of Load Expectation (LOLE),
% Loss of Energy Expectation (LOEE), Expected Interruption Rate (EIR),
% Wind Energy Share, and Expected Energy Not Served (EENS) for both
% the system without wind and with wind.
%
% Inputs:
%   load_model           - A string indicating the load modeling approach:
%                          'comulative' for cumulative load duration curve based calculations.
%                          'exact' for exact margin state calculations.
%   SCOPT_trimmed        - System Capacity-Output Probability Table, potentially trimmed.
%                          Used for calculations without wind integration.
%                          Column 1: Capacity, Column 2: Probability,
%                          Column 3 (if exists): Upward rate, Column 4 (if exists): Downward rate.
%   Annual_peak          - The annual peak load, used for trimming SCOPT_trimmed in 'comulative' model.
%   CONVOLVED_SCOPT_merged - Convolved System Capacity-Output Probability Table (with wind).
%                            Column 1: Capacity, Column 2: Probability.
%                            Column 3 (if exists): Upward rate, Column 4 (if exists): Downward rate.
%   SCOPT                - System Capacity-Output Probability Table (before trimming, if applicable).
%                          Used in 'exact' model.
%                          Column 1: Capacity, Column 2: Probability,
%                          Column 3: Upward rate, Column 4: Downward rate.
%
% Outputs:
%   LOLP                 - Loss of Load Probability.
%   LOLE                 - Loss of Load Expectation (in hours/year or days/year depending on model).
%   LOEE                 - Loss of Energy Expectation (in MWh/year).
%   EIR                  - Expected Interruption Rate.
%   Wind_Eshare          - Wind Energy Share (in %).
%   EENS_0               - Expected Energy Not Served (without wind, or total energy).
%   EENS_W               - Expected Energy Not Served (with wind).

% Default value for e (exposure factor for LOLE conversion), assuming 365 days in a year if not defined.
% This is inferred from the LOLE calculation in the 'exact' case.
e = 1; % Placeholder, as it's used as 365/e and not e itself.
       % The original code seems to imply e=1 for simplicity if 365 is used directly.
       % If 'e' is intended as a specific exposure factor, it should be passed as an input.

if strcmp(load_model, 'comulative')
    %% Cumulative Load Model Approach
    % This approach typically uses the Load Duration Curve (LDC) for reliability calculations.
    % RTS 24 load parameters cimmulative load model (hourly)
    % Weekly peak load percentages (first 52 weeks)
    weekly_peaks = [ ...
        86.2, 90, 87.8, 83.4, 88, 84.1, 83.2, 80.6, 74, 73.7, 71.5, 72.7, ...
        70.4, 75, 72.1, 80, 75.4, 83.7, 87, 88, 85.6, 81.1, 90, 88.7, 89.6, ...
        86.1, 75.5, 81.6, 80.1, 88, 72.2, 77.6, 80, 72.9, 72.6, 70.5, 78, ...
        69.5, 72.4, 72.4, 74.3, 74.4, 80, 88.1, 88.5, 90.9, 94, 89, 94.2, ...
        97, 100, 95.2];
    
    % Define week categories
    winter_weeks = [1:8, 44:52];
    summer_weeks = 18:30;
    spring_fall_weeks = [9:17, 31:43];
    
    % Hourly profiles for each type (weekday and weekend)
    profiles.winter.weekday = [67 63 60 59 59 60 74 86 95 96 96 95 ...
                               95 95 93 94 99 100 100 96 91 83 73 63];
    profiles.winter.weekend = [78 72 68 66 64 65 66 70 80 88 90 91 ...
                               90 88 87 87 91 100 99 97 94 92 87 81];
    
    profiles.summer.weekday = [64 60 58 56 56 58 64 76 87 95 99 100 ...
                               99 100 100 97 96 96 93 92 92 93 87 72];
    profiles.summer.weekend = [74 70 66 65 64 62 62 66 81 86 91 93 ...
                               93 92 91 91 92 94 95 95 100 93 88 80];
    
    profiles.spring_fall.weekday = [63 62 60 58 59 65 72 85 95 99 100 99 ...
                                     93 92 90 88 90 92 96 98 96 90 80 70];
    profiles.spring_fall.weekend = [75 73 69 66 65 65 68 74 83 89 92 94 ...
                                     91 90 90 86 85 88 92 100 97 95 90 85];
    
    % Preallocate hourly load vector
    hourly_load = zeros(52 * 7 * 24, 1);  % 52 weeks * 7 days * 24 hours
    hour_idx = 1;
    
    for w = 1:52
        % Determine season
        if ismember(w, winter_weeks)
            profile = profiles.winter;
        elseif ismember(w, summer_weeks)
            profile = profiles.summer;
        else
            profile = profiles.spring_fall;
        end
    
        % Get peak for this week
        peak = weekly_peaks(w);
    
        for d = 1:7
            if d <= 5
                daily_profile = profile.weekday;
            else
                daily_profile = profile.weekend;
            end
    
            % Scale profile by weekly peak and add to total load
            hourly_load(hour_idx:hour_idx+23) = peak * daily_profile / 100;
            hour_idx = hour_idx + 24;
        end
    end
    
    % Add day 7 of week 52 (repeat Sunday)
    w = 52;
    peak = weekly_peaks(w);
    
    if ismember(w, winter_weeks)
        profile = profiles.winter;
    elseif ismember(w, summer_weeks)
        profile = profiles.summer;
    else
        profile = profiles.spring_fall;
    end
    
    % Sunday (day 7) is a weekend
    hourly_load(hour_idx:hour_idx+23) = peak * profile.weekend / 100;
    
    % Sort descending to generate load duration curve
    sorted_load = sort(hourly_load, 'descend');
    scaled_load=sorted_load*Annual_peak/100;
    
    % Total hours in the year
    total_hours = 365 * 24;
    total_hours = linspace(0, 1, total_hours);

    LOPT=[scaled_load,total_hours.'];
    
    % % Plot with percentage x-axis
    % figure;
    % plot(total_hours, scaled_load, 'LineWidth', 1.5);
    % title('Load Duration Curve (Hourly, Percent Time)');
    % xlabel('Percentage of Time [%]');
    % ylabel('Load [% of Annual Peak]');
    % grid on;
    
    % Extract Load Duration Curve (LDC) data.
    % x_vals represents the duration (time in hours/year) that the load equals or exceeds y_vals.
    x_vals = LOPT(:,2) * 8760; % Time percentage (normalized to 1) * total hours in a year (8760)
    y_vals = LOPT(:,1);      % Scaled load (capacity)

    % Trim the SCOPT_trimmed table: only consider capacities less than the Annual_peak.
    % This ensures that only relevant capacity states for reliability assessment are used.
    SCOPT_trimmed = SCOPT_trimmed(SCOPT_trimmed(:,1) < Annual_peak, :);

    % Initialize vectors to store Loss of Energy (LE), Loss of Probability (LP),
    % and Loss of Load Expectation (LL) for each capacity state.
    LE = zeros(size(SCOPT_trimmed, 1), 1);
    LP = zeros(size(SCOPT_trimmed, 1), 1);
    LL = zeros(size(SCOPT_trimmed, 1), 1);
    
    % Calculate LE, LP, and LL for each capacity state in SCOPT_trimmed.
    % These are contributions to the overall system reliability indices.
    for i = 1:length(LE)
        threshold = SCOPT_trimmed(i,1); % Current available capacity as a threshold
    
        % Find time durations (x_vals) where the load (y_vals) is greater than or equal to the current capacity threshold.
        valid_idx = find(y_vals >= threshold);
    
        % Ensure there are at least two points for trapezoidal integration.
        if numel(valid_idx) >= 2
            x_sub = x_vals(valid_idx);        % Time durations corresponding to load exceeding threshold
            y_sub = y_vals(valid_idx) - threshold; % Load shortage for each duration
            
            % LP (Loss of Probability) for this state:
            % (Maximum duration of load exceeding threshold / total hours in year) * Probability of this capacity state.
            LP(i) = (max(x_sub) / 8760) * SCOPT_trimmed(i,2);
            
            % LL (Loss of Load Expectation) for this state:
            % (Maximum duration of load exceeding threshold) * Probability of this capacity state.
            LL(i) = (max(x_sub)) * SCOPT_trimmed(i,2);
            
            % LE (Loss of Energy Expectation) for this state:
            % Integral of (load shortage over time) * Probability of this capacity state.
            % `trapz` computes the area under the curve using the trapezoidal rule.
            LE(i) = trapz(x_sub, y_sub) * SCOPT_trimmed(i,2);
        else
            % If no load exceeds the threshold, or insufficient data points, contributions are zero.
            LE(i) = 0;
            LL(i) = 0; % Explicitly set to 0 as it's used in sum(LL)
            LP(i) = 0; % Explicitly set to 0 as it's used in sum(LP)
        end
    end

    % --- Reliability Indices (System without Wind) ---
    % 1. Loss of Load Probability (LOLP): Sum of LP contributions across all states.
    LOLP = sum(LP);
    % 2. Loss of Load Expectation (LOLE): Sum of LL contributions across all states.
    LOLE = sum(LL);
    % 3. Loss of Energy Expectation (LOEE): Sum of LE contributions.
    % Note: The original code applies a negative sign here.
    LOEE = -sum(LE);
    
    % Calculate Expected Total Energy (E_tot) from the Load Duration Curve.
    % This involves calculating the expected value of the load using its PMF.
    x_vals_load = LOPT(:,1);     % Values of the random variable (Load levels)
    cdf_vals_load = LOPT(:,2);   % Cumulative probabilities (time percentages, increasing)

    % Compute probability mass (PMF) from CDF differences.
    % The first entry is assumed to be the probability of the first load level.
    pmf_vals_load = [cdf_vals_load(1); diff(cdf_vals_load)];
    
    % Compute expected total energy by summing (load_level * probability_mass * total_hours_in_year).
    E_tot = sum(x_vals_load .* pmf_vals_load) * 8760;
    
    % Expected Interruption Rate (EIR): A measure of the expected number of interruptions per unit of energy supplied.
    % Formula: 1 + (LOEE / E_tot)
    EIR = 1 + LOEE / E_tot;

    % Initial EENS (Expected Energy Not Served) without any generation, effectively the total load energy.
    EENS_0 = E_tot;

    % --- Reliability Indices (System with Wind - Convolved System) ---
    % Initialize vectors for LE, LL, LP for the convolved system (with wind).
    LE_W = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);
    LL_W = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);
    LP_W = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);
    
    % Initialize vectors for "Surplus" Energy (LE_S, LL_S, LP_S).
    % These capture instances where generation exceeds load.
    LE_S = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);
    LL_S = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);
    LP_S = zeros(size(CONVOLVED_SCOPT_merged, 1), 1);

    % Perform calculations for the convolved system (WPP + MTDC).
    for i = 1:length(LE_W)
        threshold = CONVOLVED_SCOPT_merged(i,1); % Current available capacity of the convolved system
    
        % --- Calculate Loss of Energy/Load/Probability (for deficit) ---
        % Find time durations where load (y_vals) is greater than or equal to the convolved capacity threshold.
        valid_idx = find(y_vals >= threshold);
    
        if numel(valid_idx) >= 2
            x_sub = x_vals(valid_idx);
            y_sub = y_vals(valid_idx) - threshold; % Load shortage
            
            LE_W(i) = trapz(x_sub, y_sub) * CONVOLVED_SCOPT_merged(i,2); % Energy not served contribution
            LP_W(i) = (max(x_sub) / 8760) * CONVOLVED_SCOPT_merged(i,2); % Loss probability contribution
            LL_W(i) = (max(x_sub)) * CONVOLVED_SCOPT_merged(i,2); % Loss of load expectation contribution
        else
            % If no load exceeds the threshold, or insufficient data points, contributions are zero.
            LE_W(i) = 0;
            LL_W(i) = 0;
            LP_W(i) = 0;
        end

        % --- Calculate Surplus Energy/Load/Probability ---
        % Find time durations where load (y_vals) is less than the convolved capacity threshold.
        % This represents instances where generation capacity exceeds demand.
        valid_idx_surplus = find(y_vals < threshold);

        if numel(valid_idx_surplus) >= 2
            x_sub_s = x_vals(valid_idx_surplus);
            y_sub_s = -y_vals(valid_idx_surplus) + threshold; % Surplus energy (positive value)
            
            LE_S(i) = trapz(x_sub_s, y_sub_s) * CONVOLVED_SCOPT_merged(i,2); % Surplus energy contribution
            LP_S(i) = (max(x_sub_s) / 8760) * CONVOLVED_SCOPT_merged(i,2); % Surplus probability contribution
            LL_S(i) = (max(x_sub_s)) * CONVOLVED_SCOPT_merged(i,2); % Surplus load expectation contribution
        else
            LE_S(i) = 0;
            LL_S(i) = 0;
            LP_S(i) = 0;
        end
    end
    
    % EENS_W (Expected Energy Not Served with Wind):
    % Calculated as the negative sum of LE_W contributions.
    % Note: The negative sign indicates it's an energy deficit.
    EENS_W = -sum(LE_W); 
    
    % Wind_Eshare: Percentage share of wind energy in meeting demand.
    % Formula: 100 * (EENS_0 + EENS_W) / (EENS_0 * EIR)
    % This formula appears specific to the context and combines EENS with and without wind, and EIR.
    Wind_Eshare = 100 * (EENS_0 + EENS_W) / (EENS_0 * EIR);

elseif strcmp(load_model, 'exact')
    %% Exact Load Model Approach
    % This approach calculates margins for all combinations of supply and demand states.
    % RTS 24 load parameters exact load model
    % Load model RTS 24 
    %   'Load X (MW)   Prob. p(X)       λ+ (occ/day)      λ- (occ/day)      freq. (occ/day)' 
    LOPT=[  2687	    0.016438		    0	                2	            0.032877
            2454	    0.112329	        0	                2           	0.224658;
            2188	    0.147945	        0	                2           	0.032877;
            1953	    0.158904	        0	                2           	0.317808;
            1593	    0.064384	        0	                2           	0.032877;
            1485	    0.500000	        2	                0           	1];
    % --- Calculations for System without Wind ---
    % Number of System Capacity states and Load states.
    nSC = size(SCOPT, 1); % Using SCOPT (untrimmed) for exact model
    nLO = size(LOPT, 1);

    % Initialize marginTable to store combined margin states.
    % Columns: [Margin, Probability, Lambda_plus, Lambda_minus, Frequency]
    marginTable = zeros(nSC * nLO, 5);
    
    index = 0; % Index for populating marginTable
    % Iterate through all combinations of SCOPT and LOPT states.
    for i = 1:nSC % Loop through SCOPT states
        for j = 1:nLO % Loop through LOPT states
            index = index + 1;
            % Margin: Available Capacity - Load
            margin = SCOPT(i, 1) - LOPT(j, 1);
            % Combined Probability: Product of individual probabilities (assuming independence).
            probability = SCOPT(i, 2) * LOPT(j, 2);
            % Combined Upward Transition Rate: Sum of SCOPT's upward rate and LOPT's downward rate.
            lambda_plus = SCOPT(i, 3) + LOPT(j, 4);
            % Combined Downward Transition Rate: Sum of SCOPT's downward rate and LOPT's upward rate.
            lambda_minus = SCOPT(i, 4) + LOPT(j, 3);
            % Combined Frequency: Probability * (Lambda_plus + Lambda_minus)
            frequency = probability * (lambda_plus + lambda_minus);
            
            marginTable(index, :) = [margin, probability, lambda_plus, lambda_minus, frequency];
        end
    end

    % Combine identical margin states in marginTable to form MOPT (Margin-Output Probability Table).
    uniqueMargins = unique(marginTable(:, 1));
    MOPT = zeros(length(uniqueMargins), 5);
    
    for k = 1:length(uniqueMargins)
        margin = uniqueMargins(k);
        indices = find(marginTable(:, 1) == margin); % Find all entries with this unique margin
        
        % Aggregate probabilities and frequencies for the current unique margin.
        aggregatedProb = sum(marginTable(indices, 2));
        aggregatedFreq = sum(marginTable(indices, 5));
        
        % Calculate aggregated lambda_plus and lambda_minus using Equation 3.26 (weighted average).
        % This is a probability-weighted average of the individual lambda_plus/minus for states within this margin.
        lambda_plus_aggr = sum(marginTable(indices, 2) .* marginTable(indices, 3));
        if aggregatedProb ~= 0
            lambda_plus_aggr = lambda_plus_aggr / aggregatedProb;
        else
            lambda_plus_aggr = 0; % Avoid division by zero
        end

        lambda_minus_aggr = sum(marginTable(indices, 2) .* marginTable(indices, 4));
        if aggregatedProb ~= 0
            lambda_minus_aggr = lambda_minus_aggr / aggregatedProb;
        else
            lambda_minus_aggr = 0; % Avoid division by zero
        end
        
        MOPT(k, :) = [margin, aggregatedProb, lambda_plus_aggr, lambda_minus_aggr, aggregatedFreq];
    end
    
    % --- Reliability Indices (System without Wind) ---
    % Identify states where margin is negative (Load > Capacity).
    negativeMargins = MOPT(:, 1) < 0; 
    
    % 1. Loss of Load Probability (LOLP): Sum of probabilities for all negative margin states.
    LOLP = sum(MOPT(negativeMargins, 2)); 
    
    % 2. Loss of Load Expectation (LOLE):
    % The original code calculates LOLE as (probability for negative margin) * 365 / e.
    % If `e` is an exposure factor, this would convert probability to days/year of loss.
    % The subsequent `sum(LOLE)` suggests LOLE is summed over all negative margin states.
    LOLE_individual_states = MOPT(negativeMargins, 2) * 365 / e; 
    LOLE = sum(LOLE_individual_states); % Sum across all negative margins

    % 3. Loss of Energy Expectation (LOEE):
    % Calculated as (negative margin) * (probability for negative margin) * 365 * 24 (hours in a year).
    % This represents the total energy not supplied due to capacity shortage.
    LOEE_individual_states = MOPT(negativeMargins, 1) .* MOPT(negativeMargins, 2) * 365 * 24; 
    LOEE = sum(LOEE_individual_states); % Sum across all negative margins
    
    % Store LOEE (if EENS(x,y) was part of a larger context, it would be assigned here)
    % EENS(x,y) = LOEE; % This line is commented out as 'x' and 'y' are not defined in this function scope.
    
    % E_tot: Total expected energy demanded by the load over a year.
    E_tot = sum(LOPT(:,1) .* LOPT(:,2)) * 365 * 24;
    
    % Expected Interruption Rate (EIR):
    EIR = 1 + LOEE / E_tot;

    % EENS_0: Expected Energy Not Served (without wind), here defined as E_tot.
    EENS_0 = E_tot;   
    
    % --- Calculations for System with Wind (Convolved System) ---
    % Similar margin analysis, but using the CONVOLVED_SCOPT_merged data.
    nSC_conv = size(CONVOLVED_SCOPT_merged, 1);
    nLO_conv = size(LOPT, 1);
    marginTable_W = zeros(nSC_conv * nLO_conv, 5);
    
    index_W = 0;
    for i = 1:nSC_conv
        for j = 1:nLO_conv
            index_W = index_W + 1;
            margin_W = CONVOLVED_SCOPT_merged(i, 1) - LOPT(j, 1);
            probability_W = CONVOLVED_SCOPT_merged(i, 2) * LOPT(j, 2);
            % Note: Assuming LOPT has 4 columns (capacity, prob, lambda_plus, lambda_minus) for these rates.
            lambda_plus_W = CONVOLVED_SCOPT_merged(i, 3) + LOPT(j, 4); 
            lambda_minus_W = CONVOLVED_SCOPT_merged(i, 4) + LOPT(j, 3);
            frequency_W = probability_W * (lambda_plus_W + lambda_minus_W);
            
            marginTable_W(index_W, :) = [margin_W, probability_W, lambda_plus_W, lambda_minus_W, frequency_W];
        end
    end
    
    % Combine identical margin states for the wind-integrated system.
    uniqueMargins_W = unique(marginTable_W(:, 1));
    MOPT_W = zeros(length(uniqueMargins_W), 5);
    
    for k = 1:length(uniqueMargins_W)
        margin_W = uniqueMargins_W(k);
        indices_W = find(marginTable_W(:, 1) == margin_W);
        aggregatedProb_W = sum(marginTable_W(indices_W, 2));
        aggregatedFreq_W = sum(marginTable_W(indices_W, 5));
        
        lambda_plus_aggr_W = sum(marginTable_W(indices_W, 2) .* marginTable_W(indices_W, 3));
        if aggregatedProb_W ~= 0
            lambda_plus_aggr_W = lambda_plus_aggr_W / aggregatedProb_W;
        else
            lambda_plus_aggr_W = 0;
        end

        lambda_minus_aggr_W = sum(marginTable_W(indices_W, 2) .* marginTable_W(indices_W, 4));
        if aggregatedProb_W ~= 0
            lambda_minus_aggr_W = lambda_minus_aggr_W / aggregatedProb_W;
        else
            lambda_minus_aggr_W = 0;
        end
        
        MOPT_W(k, :) = [margin_W, aggregatedProb_W, lambda_plus_aggr_W, lambda_minus_aggr_W, aggregatedFreq_W];
    end
    
    % RE-CALCULATE negativeMargins SPECIFICALLY FOR MOPT_W
    % This is the fix for the reported error.
    negativeMargins_W = MOPT_W(:, 1) < 0;

    % EENS_W (Expected Energy Not Served with Wind):
    % Calculate EENS based on negative margins in MOPT_W.
    EENS_W = sum(MOPT_W(negativeMargins_W, 1) .* MOPT_W(negativeMargins_W, 2)) * 365 * 24;
    
    % Wind_Eshare: Percentage share of wind energy.
    Wind_Eshare = 100 * (EENS_0 + EENS_W) / (EENS_0 * EIR);

else
    % Error handling for invalid load_model input.
    error('Invalid load_model specified. Choose either ''comulative'' or ''exact''.');
end

end