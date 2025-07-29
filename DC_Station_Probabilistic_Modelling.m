function [SEDCTCPCOPT, REDCTCPCOPT, DCTCP_TM_Reduced, DCTCPCOPT] = DC_Station_Probabilistic_Modelling( ...
    MTTR_ACB, MTTR_CT, MTTR_CPD, MTTR_ACF, MTTR_CR, MTTR_DCL, ...
    lambda_ACB, lambda_CT, lambda_CPD, lambda_ACF, lambda_CR, ...
    lambda_arm_optimum, lambda_con, A_mmc_optimum, ...
    lambda_DCB, mu_DCB, lambda_DCL, cable_length, ...
    WPP_CAP, MTDC_CAP, nSEDCTCPs, nREDCTCPs)
%DC_Station_Probabilistic_Modelling Calculates probabilistic metrics for a DC transmission station.
%
%   [SEDCTCPCOPT, REDCTCPCOPT, DCTCP_TM] = DC_Station_Probabilistic_Modelling(...)
%
%   This function models the reliability of a DC transmission station by
%   constructing its state transition matrix, calculating limiting state
%   probabilities, and deriving capacity outage probability tables for
%   both the sending and receiving ends.
%
%   Inputs:
%   - MTTR_ACB (double): Mean Time To Repair for AC Breaker (hours).
%   - MTTR_CT (double): Mean Time To Repair for Converter Transformer (hours).
%   - MTTR_CPD (double): Mean Time To Repair for Converter Protection Device (hours).
%   - MTTR_ACF (double): Mean Time To Repair for AC Filter (hours).
%   - MTTR_CR (double): Mean Time To Repair for Converter Reactor (hours).
%   - MTTR_DCL (double): Mean Time To Repair for DC Cable Link (hours/100km).
%   - lambda_ACB (double): Failure rate of AC Breaker (occ/y).
%   - lambda_CT (double): Failure rate of Converter Transformer (occ/y).
%   - lambda_CPD (double): Failure rate of Converter Protection Device (occ/y).
%   - lambda_ACF (double): Failure rate of AC Filter (occ/y).
%   - lambda_CR (double): Failure rate of Converter Reactor (occ/y).
%   - lambda_arm_optimum (double): Optimal failure rate of a single MMC arm (from MMC_Probabilistic_Modelling) (occ/y).
%   - lambda_con (double): Failure rate at the converter level (occ/y).
%   - A_mmc_optimum (double): Optimal availability of the overall MMC (from MMC_Probabilistic_Modelling).
%   - lambda_DCB (double): Failure rate of DC Breaker (occ/y).
%   - mu_DCB (double): Repair rate of DC Breaker (occ/y).
%   - lambda_DCL (double): Failure rate of DC Cable Link (occ/y per 100km).
%   - cable_length (double): Length of the DC cable in km.
%   - WPP_CAP (double): Wind Power Plant Capacity in MW.
%   - MTDC_CAP (double): Multi-Terminal DC Capacity in MW (used for rounding/scaling).
%   - nSEDCTCPs (double): Number of Sending End DC Transmission Capacity Parts.
%   - nREDCTCPs (double): Number of Receiving End DC Transmission Capacity Parts.
%
%   Outputs:
%   - SEDCTCPCOPT (matrix): Sending End DC Transmission Capacity Outage Probability Table.
%     Columns: [Capacity (MW), Probability, Upward Transition Rate,
%     Downward Transition Rate, Frequency].
%   - REDCTCPCOPT (matrix): Receiving End DC Transmission Capacity Outage Probability Table.
%     (Structure identical to SEDCTCPCOPT).
%   - DCTCP_TM_Reduced (matrix): 3-States DC Transmission Capacity Probability Transition Matrix (occ/day).
%   - DCTCPCOPT (matrix): % Generic 3-states DC Transmission Capacity Probability probability table.

    %% Repair Rate Definitions
    % Convert Mean Time To Repair (MTTR) from hours to occurrences per year.
    % Note: mu_DCB is given directly in (occ/y).
    mu_ACB = 8760 / MTTR_ACB;   % AC Breaker repair rate (occ/y)
    mu_CT = 8760 / MTTR_CT;     % Converter Transformer repair rate (occ/y)
    mu_CPD = 8760 / MTTR_CPD;   % Converter Protection Device repair rate (occ/y)
    mu_ACF = 8760 / MTTR_ACF;   % AC Filter repair rate (occ/y)
    mu_CR = 8760 / MTTR_CR;     % Converter Reactor repair rate (occ/y)
    mu_DCL = 8760 / MTTR_DCL;   % DC Cable Link repair rate (occ/y per 100km)
    
    %% Aggregated Failure and Repair Rates for Sub-systems
    % These rates represent the effective failure/repair rates for lumped
    % sub-systems (station components, MMC, and DC cable).
    
    % Sub-system 'a': Aggregated AC-side station components
    lambda_sa = lambda_ACB + lambda_CT + lambda_CPD + lambda_ACF + lambda_CR;
    mu_sa = lambda_sa / (lambda_ACB/mu_ACB + lambda_CT/mu_CT + lambda_CPD/mu_CPD + lambda_ACF/mu_ACF + lambda_CR/mu_CR);
    
    % Sub-system 'b': Aggregated MMC (Modular Multilevel Converter)
    % lambda_arm_optimum and A_mmc_optimum are assumed to come from the MMC_Probabilistic_Modelling.
    lambda_sb = 6 * lambda_arm_optimum + lambda_con; % Assuming 6 arms in the MMC
    mu_sb = (lambda_sb / (1 - A_mmc_optimum)) - lambda_sb; % Derived from availability definition A = mu / (lambda + mu)
    
    % Sub-system 'c': Aggregated DC-side (DC Breaker and DC Cable Link)
    % Cable length is divided by 2 as the failure rate given is per 100km for the entire link,
    % but the model seems to imply two segments for a full length.
    lambda_sc = lambda_DCB + lambda_DCL * cable_length / (2 * 100);
    mu_sc = lambda_sc / (lambda_DCB/mu_DCB + (lambda_DCL * cable_length / (2 * 100)) / (mu_DCL * cable_length / (2 * 100)));
    
    %% DC Transmission Capacity Probability Transition Matrix (DCTCP_TM)
    % This section builds the full transition matrix for the 13 states
    % of the DC Transmission Capacity system. The states represent various
    % combinations of operational and failed sub-systems.
    % The matrix is initially populated with rates in (occ/y) and then converted to (occ/day).
    
    % Initialize the transition matrix (13x13 states)
    DCTCP_TM = zeros(13, 13);
    
    % Fill the transition matrix based on specified state transitions
    % The states (1 to 13) implicitly represent different operational/outage
    % conditions of the aggregated sub-systems 'a', 'b', and 'c'.
    % E.g., a transition from state 1 to state 2 means a failure of sub-system 'b'.
    
    % Transitions FROM State 1 (Fully Operational)
    DCTCP_TM(1, 2) = 2 * lambda_sb; % From 1 (0,0,0) to state with 1 'b' failure
    DCTCP_TM(1, 3) = 2 * lambda_sc; % From 1 (0,0,0) to state with 1 'c' failure
    DCTCP_TM(1, 5) = lambda_sa;     % From 1 (0,0,0) to state with 1 'a' failure
    
    % Transitions FROM State 2 (e.g., 1 'b' failure)
    DCTCP_TM(2, 1) = mu_sb;         % Repair of 'b'
    DCTCP_TM(2, 4) = 2 * lambda_sc; % Failure of 2 'c's
    DCTCP_TM(2, 6) = lambda_sb;     % Another 'b' failure
    DCTCP_TM(2, 7) = lambda_sa;     % 'a' failure
    DCTCP_TM(2, 10) = 2 * lambda_sc; % Failure of 2 'c's
    
    % Transitions FROM State 3 (e.g., 1 'c' failure)
    DCTCP_TM(3, 1) = mu_sc;         % Repair of 'c'
    DCTCP_TM(3, 4) = 2 * lambda_sb; % Failure of 2 'b's
    DCTCP_TM(3, 8) = lambda_sc;     % Another 'c' failure
    DCTCP_TM(3, 9) = lambda_sa;     % 'a' failure
    DCTCP_TM(3, 10) = 2 * lambda_sb; % Another 'b' failure
    
    % Transitions FROM State 4 (e.g., 1 'b' failure and 1 'c' failure)
    DCTCP_TM(4, 2) = mu_sc;         % Repair of 'c'
    DCTCP_TM(4, 3) = mu_sb;         % Repair of 'b'
    DCTCP_TM(4, 11) = lambda_sa;    % 'a' failure
    DCTCP_TM(4, 12) = lambda_sc;    % Another 'c' failure
    DCTCP_TM(4, 13) = lambda_sb;    % Another 'b' failure
    
    % Transitions FROM State 5 (e.g., 1 'a' failure)
    DCTCP_TM(5, 1) = mu_sa;         % Repair of 'a'
    
    % Transitions FROM State 6 (e.g., 2 'b' failures)
    DCTCP_TM(6, 2) = 2 * mu_sb;     % Repair of 2 'b's
    
    % Transitions FROM State 7 (e.g., 1 'a' failure and 1 'b' failure)
    DCTCP_TM(7, 2) = mu_sa;         % Repair of 'a'
    
    % Transitions FROM State 8 (e.g., 2 'c' failures)
    DCTCP_TM(8, 3) = 2 * mu_sc;     % Repair of 2 'c's
    
    % Transitions FROM State 9 (e.g., 1 'a' failure and 1 'c' failure)
    DCTCP_TM(9, 3) = mu_sa;         % Repair of 'a'
    
    % Transitions FROM State 10 (e.g., 1 'b' failure and 2 'c' failures)
    DCTCP_TM(10, 2) = mu_sc;        % Repair of 'c'
    DCTCP_TM(10, 3) = mu_sb;        % Repair of 'b'
    
    % Transitions FROM State 11 (e.g., 1 'a' failure, 1 'b' failure, and 1 'c' failure)
    DCTCP_TM(11, 4) = mu_sa;        % Repair of 'a'
    
    % Transitions FROM State 12 (e.g., 1 'b' failure, 1 'c' failure, and another 'c' failure)
    DCTCP_TM(12, 4) = 2 * mu_sc;    % Repair of 2 'c's
    
    % Transitions FROM State 13 (e.g., 1 'b' failure, 1 'c' failure, and another 'b' failure)
    DCTCP_TM(13, 4) = 2 * mu_sb;    % Repair of 2 'b's
    
    % Convert transition matrix from (occ/y) to (occ/day)
    DCTCP_TM = DCTCP_TM / 365;

    %% Limiting State Probabilities of DCTCP
    % Calculates the steady-state probabilities of the DC transmission system
    % being in each of its 13 defined states. This is done by finding the
    % eigenvector corresponding to an eigenvalue of 1 for the transition matrix.
    
    % Define a time step (1 year). This Delta_t is conceptual for building
    % the State Transition Probability (STP) matrix from rates.
    Delta_t = 1; % years
    
    % Create the State Transition Probability (STP) matrix
    % STP(i,j) represents the probability of transitioning from state i to state j in Delta_t.
    STP_DCTCP_TM = DCTCP_TM * Delta_t;
    
    % Adjust the diagonal entries of the STP matrix. The sum of probabilities
    % from any state must be 1. The diagonal represents the probability of
    % remaining in the same state.
    for i = 1:size(STP_DCTCP_TM, 1)
        STP_DCTCP_TM(i, i) = 1 - sum(STP_DCTCP_TM(i, :)) + STP_DCTCP_TM(i, i); % Add back the original diagonal to subtract full row sum
    end
    
    % Handle potential issues with sum for diagonal (should be 1 - sum(off-diagonals))
    % The correct calculation for continuous-time Markov chains to get the
    % generator matrix Q from transition rates is to set diagonal elements Q(i,i) = -sum(Q(i,j) for j!=i).
    % Then the steady-state probabilities alpha * Q = 0.
    % However, the previous line `STP_DCTCP_TM(i,i) = 1 - sum(STP_DCTCP_TM(i,:))`
    % effectively creates a row-stochastic matrix if DCTCP_TM was non-negative.
    % Given the initial `DCTCP_TM = zeros(13,13)` and then positive rates,
    % the `STP_DCTCP_TM = DCTCP_TM * Delta_t;` followed by `STP_DCTCP_TM(i,i) = 1 - sum(STP_DCTCP_TM(i,:));`
    % is a common way to form the one-step transition probability matrix P
    % from a rate matrix (Q=P-I, where I is identity matrix and P is probability matrix).
    % Let's re-verify the diagonal adjustment for clarity:
    for i = 1:size(STP_DCTCP_TM, 1)
        STP_DCTCP_TM(i, i) = 1 - sum(DCTCP_TM(i, :)); % Sum of off-diagonal rates for state i
    end

    % Calculate the limiting state probabilities (alpha_DCTCP) by finding the
    % eigenvector corresponding to the eigenvalue of 1 of the transpose of the STP matrix.
    % This is based on the equation pi * P = pi, where pi is the vector of steady-state probabilities.
    [V, D] = eig(STP_DCTCP_TM.');
    
    % Find the index of the eigenvalue that is closest to 1
    [~, idx] = min(abs(diag(D) - 1));
    
    % Extract the corresponding eigenvector
    alpha_DCTCP = V(:, idx);
    
    % Normalize the eigenvector so that its elements sum to 1,
    % ensuring it represents probabilities.
    alpha_DCTCP = alpha_DCTCP / sum(alpha_DCTCP);
    
    %% Building DC Transmission Capacity Probability Tables (DCTCPCOPT)
    % This section translates the 13 detailed system states into three
    % aggregated capacity states (100%, 50%, 0%) for both sending and
    % receiving ends, populating the SEDCTCPCOPT and REDCTCPCOPT tables.
    
    % Calculate the single block capacity for sending and receiving ends.
    % This rounds the WPP_CAP to the nearest MTDC_CAP multiple and divides
    % by the number of capacity parts.
    SDCT_cap_base = round(WPP_CAP / MTDC_CAP) * MTDC_CAP / nSEDCTCPs;
    RDCT_cap_base = round(WPP_CAP / MTDC_CAP) * MTDC_CAP / nREDCTCPs;
    
    % Define the actual capacity states (100%, 50%, 0%) for a single block
    SDCT_cap_states = SDCT_cap_base * [1, 0.5, 0];
    RDCT_cap_states = RDCT_cap_base * [1, 0.5, 0];
    
    % Initialize the Sending End DCTCPCOPT table.
    % Columns: [Capacity (MW), Probability, Upward Rate, Downward Rate, Frequency]
    SEDCTCPCOPT = zeros(length(SDCT_cap_states), 5);
    SEDCTCPCOPT(:, 1) = SDCT_cap_states'; % Set capacity column
    
    % Assign probabilities to the aggregated capacity states based on alpha_DCTCP.
    % Assuming a mapping of the 13 states to the 3 aggregated states:
    % State 1: 100% capacity
    % States 2, 3, 4: 50% capacity
    % States 5 to 13: 0% capacity
    SEDCTCPCOPT(1, 2) = alpha_DCTCP(1);             % Probability for 100% capacity
    SEDCTCPCOPT(2, 2) = sum(alpha_DCTCP(2:4));      % Probability for 50% capacity
    SEDCTCPCOPT(3, 2) = sum(alpha_DCTCP(5:13));     % Probability for 0% capacity
    
    % Extract probabilities for easier use in subsequent calculations
    p_k = SEDCTCPCOPT(:, 2);
    
    % Initialize arrays for + and - transition rates for each aggregated state
    lambda_plus = zeros(3, 1);  % Upward transition rates (to higher capacity)
    lambda_minus = zeros(3, 1); % Downward transition rates (to lower capacity)
    
    % Calculate aggregated upward (lambda_plus) and downward (lambda_minus)
    % transition rates between the three capacity states. This involves summing
    % the appropriate elements from the original DCTCP_TM.
    
    % State 1: 100% Capacity (Original state 1)
    i_range_100 = 1;
    % Transitions from 100% to 50% (States 2,3,4) or 0% (States 5-13) are downward.
    lambda_minus(1) = sum(sum(DCTCP_TM(i_range_100, [2:13]), 2));
    lambda_plus(1) = 0; % Cannot transition to higher capacity from 100%
    
    % State 2: 50% Capacity (Original states 2, 3, 4)
    i_range_50 = 2:4;
    % Transitions from 50% to 100% (State 1) are upward.
    lambda_plus(2) = sum(sum(DCTCP_TM(i_range_50, 1), 2));
    % Transitions from 50% to 0% (States 5-13) are downward.
    lambda_minus(2) = sum(sum(DCTCP_TM(i_range_50, [5:13]), 2));
    
    % State 3: 0% Capacity (Original states 5 to 13)
    i_range_0 = 5:13;
    % Transitions from 0% to 100% (State 1) or 50% (States 2,3,4) are upward.
    lambda_plus(3) = sum(sum(DCTCP_TM(i_range_0, [1:4]), 2));
    lambda_minus(3) = 0; % Cannot transition to lower capacity from 0%
    
    % Calculate the frequency of each aggregated state
    % Frequency (f_k) = Probability (p_k) * (lambda_plus + lambda_minus)
    f_k = p_k .* (lambda_plus + lambda_minus);
    
    % Fill the remaining columns of the SEDCTCPCOPT table
    SEDCTCPCOPT(:, 3) = lambda_plus;  % Upward transition rates
    SEDCTCPCOPT(:, 4) = lambda_minus; % Downward transition rates
    SEDCTCPCOPT(:, 5) = f_k;          % Frequencies
    
    % The Receiving End table is assumed to have the same 3-states probabilistic
    % characteristics (probabilities, transition rates, frequencies) as the
    % Sending End, but with its own defined capacity values.
    REDCTCPCOPT = zeros(length(RDCT_cap_states), 5);
    REDCTCPCOPT(:, 1) = RDCT_cap_states'; % Set receiving end capacity column
    REDCTCPCOPT(:, 2:end) = SEDCTCPCOPT(:, 2:end); % Copy probabilities, rates, frequencies
    DCTCPCOPT = [SEDCTCPCOPT(:,1)/SDCT_cap_base,SEDCTCPCOPT(:,2)]; % Generic three states station probability table

    %% Building equivalent 3 state transition matrix in (occ/day)
        DCTCP_TM_Reduced=zeros(3,3);
        DCTCP_TM_Reduced(1,2)= sum(DCTCP_TM(1, 2:4)*alpha_DCTCP(1))/DCTCPCOPT(1,2);DCTCP_TM_Reduced(1,3)= sum(DCTCP_TM(1, 5:end)*alpha_DCTCP(1))/DCTCPCOPT(1,2);
        DCTCP_TM_Reduced(2,1)= sum(sum(alpha_DCTCP(2:4)'*DCTCP_TM(2:4, 1), 2))/DCTCPCOPT(2,2);DCTCP_TM_Reduced(2,3)= sum(sum(alpha_DCTCP(2:4)'*DCTCP_TM(2:4, 5:end), 2))/DCTCPCOPT(2,2); 
        DCTCP_TM_Reduced(3,1)= sum(sum(alpha_DCTCP(5:end)'*DCTCP_TM(5:end, 1), 2))/DCTCPCOPT(3,2);DCTCP_TM_Reduced(3,2)= sum(sum(alpha_DCTCP(5:end)'*DCTCP_TM(5:end, 2:4), 2))/DCTCPCOPT(3,2);
        

end