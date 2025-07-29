function [MTDCCOPT_merged, MTDC_TM_merged] = MTDC_Probabilistic_Modelling(nSEDCTCPs, nREDCTCPs, DCTCPCOPT, DCTCP_TM, MTDC_CAP)
% MTDC_Probabilistic_Modelling performs a probabilistic reliability assessment
% and state merging for Multi-Terminal DC (MTDC) systems. It models the
% system based on the states of its individual DC/DC Converter Transformer
% Converters (DCTCPs), calculates system capacities, state probabilities,
% and transition rates, finally merging similar states for a more compact representation.
%
% Inputs:
%   nSEDCTCPs   - Number of Sending End DCTCPs. These contribute to the sending capacity.
%   nREDCTCPs   - Number of Receiving End DCTCPs. These contribute to the receiving capacity.
%   DCTCPCOPT   - A 3x2 matrix. Row 1: N (Normal/100% capacity), Row 2: P (Partial/50% capacity),
%                 Row 3: F (Failed/0% capacity). Column 1 contains capacity, Column 2 contains probability.
%                 Example: [1, prob_N; 0.5, prob_P; 0, prob_F]
%   DCTCP_TM    - A 3x3 transition matrix for a single DCTCP. Rows and columns
%                 correspond to states N, P, F. DCTCP_TM(i, j) is the transition
%                 rate from state i to state j.
%   MTDC_CAP    - The maximum nominal capacity of the entire MTDC system (e.g., in MW).
%                 Used to scale the normalized capacities to actual values.
%
% Outputs:
%   MTDCCOPT_merged - A matrix containing the merged states' information:
%                     Column 1: Merged capacity (actual MW, rounded).
%                     Column 2: Merged probability of being in this capacity state.
%                     Column 3: Merged upward transition rate (from this state to a higher capacity state).
%                     Column 4: Merged downward transition rate (from this state to a lower capacity state).
%                     Column 5: Total departure rate from this merged state (sum of transition rates * probability).
%   MTDC_TM_merged  - A square matrix representing the merged transition rates
%                     between the merged capacity states. MTDC_TM_merged(i, j)
%                     is the transition rate from merged state i to merged state j.

        %% 1. System State Enumeration
        % Defines the possible operational states for each individual DCTCP.
        % 'N' = Normal (100% capacity), 'P' = Partial (50% capacity), 'F' = Failed (0% capacity).
        states = ['N', 'P', 'F']; 
        
        % Calculate the total number of DCTCPs in the system (sending + receiving).
        num_MTDC_ends = nSEDCTCPs + nREDCTCPs;
        % Determine the number of unique states per DCTCP (3: N, P, F).
        num_states = length(states);
        % Calculate the total number of possible combined system states.
        % This is (number of states per DCTCP) ^ (total number of DCTCPs).
        num_combinations = num_states^num_MTDC_ends;
        
        % Initialize a cell array to store all possible combinations of DCTCP states.
        % Each row represents a unique system state, and columns represent individual DCTCPs.
        ternary_table = cell(num_combinations, num_MTDC_ends);
        
        % Generate all combinations of states. This loop effectively creates a
        % base-3 representation for each combination, where each digit corresponds
        % to a DCTCP's state (0=N, 1=P, 2=F, or similar mapping).
        for i = 1:num_combinations
            temp = i - 1; % Adjust to 0-based indexing for modulo operations
            for j = num_MTDC_ends:-1:1 % Populate from rightmost DCTCP to leftmost
                % Calculate the state index for the current DCTCP (mod 3 gives 0, 1, or 2)
                state_idx = mod(temp, num_states) + 1; % Add 1 to get 1-based index for 'states' array
                % Assign the corresponding state character ('N', 'P', 'F') to the table
                ternary_table{i, j} = states(state_idx);
                % Move to the next "digit" by integer division
                temp = floor(temp / num_states);
            end
        end
        
        %% 2. Capacity Mapping for Sending and Receiving Ends
        % Calls helper functions to generate capacity mappings for the sending
        % and receiving DCTCPs based on their individual states.
        % 'sending_capacity_mapping': A map where keys are string combinations
        %                             of sending DCTCP states (e.g., 'NN', 'NP')
        %                             and values are their corresponding normalized capacities.
        % 'SDCT_cap_states_raw': All unique normalized capacity values for sending DCTCPs.
        [sending_capacity_mapping, SDCT_cap_states_raw] = MTDC_functions_helper(nSEDCTCPs);
        % Same for receiving DCTCPs.
        [receiving_capacity_mapping, RDCT_cap_states_raw] = MTDC_functions_helper(nREDCTCPs);

        % Extract unique and sorted capacity states for merging later.
        SDCT_cap_states = unique(SDCT_cap_states_raw);
        RDCT_cap_states = unique(RDCT_cap_states_raw);
        
        % Initialize a vector to store the combined system capacity for each
        % possible system state defined in 'ternary_table'.
        capacity_vector = zeros(num_combinations, 1);
        
        % Calculate the overall system capacity for each combination of DCTCP states.
        % The system capacity is the minimum of the total sending and total receiving capacities.
        for i = 1:num_combinations
            % Extract the states for sending DCTCPs (first nSEDCTCPs columns).
            sending_state_cells = ternary_table(i, 1:nSEDCTCPs);
            % Extract the states for receiving DCTCPs (remaining columns).
            receiving_state_cells = ternary_table(i, nSEDCTCPs+1:num_MTDC_ends);
            
            % Convert cell arrays of characters to a single string for map lookup.
            sending_state_str = strcat(sending_state_cells{:});
            receiving_state_str = strcat(receiving_state_cells{:});
            
            % Initialize current capacities to 0.
            current_sending_capacity = 0;
            current_receiving_capacity = 0;

            % Determine sending capacity:
            % If there are sending DCTCPs, look up their combined capacity.
            % If nSEDCTCPs is 0, assume full sending capacity (normalized to 1).
            if nSEDCTCPs > 0
                if isKey(sending_capacity_mapping, sending_state_str)
                    current_sending_capacity = sending_capacity_mapping(sending_state_str);
                end
            else
                current_sending_capacity = 1; % No sending components means no sending limitation
            end

            % Determine receiving capacity:
            % If there are receiving DCTCPs, look up their combined capacity.
            % If nREDCTCPs is 0, assume full receiving capacity (normalized to 1).
            if nREDCTCPs > 0
                if isKey(receiving_capacity_mapping, receiving_state_str)
                    current_receiving_capacity = receiving_capacity_mapping(receiving_state_str);
                end
            else
                current_receiving_capacity = 1; % No receiving components means no receiving limitation
            end
            
            % The overall system capacity is limited by the minimum of sending and receiving capacities.
            capacity_vector(i) = min(current_sending_capacity, current_receiving_capacity);
        end
        
        %% 3. Probability Calculation for Each System State
        % Initialize a vector to store the probability of each combined system state.
        probability_vector = zeros(num_combinations, 1);
        
        % Compute the joint probability for each combination of DCTCP states.
        % This is done by multiplying the individual probabilities of each DCTCP's state.
        for i = 1:num_combinations
            combination_prob = 1; % Start with a probability of 1 for the current combination
            for j = 1:num_MTDC_ends
                state = ternary_table{i, j}; % Get the state of the j-th DCTCP in the i-th combination
                % Multiply by the corresponding probability from DCTCPCOPT based on the state.
                if strcmp(state, 'N')
                    combination_prob = combination_prob * DCTCPCOPT(1,2); % Probability of Normal state
                elseif strcmp(state, 'P')
                    combination_prob = combination_prob * DCTCPCOPT(2,2); % Probability of Partial state
                elseif strcmp(state, 'F')
                    combination_prob = combination_prob * DCTCPCOPT(3,2); % Probability of Failed state
                end
            end
            probability_vector(i) = combination_prob; % Store the calculated probability
        end
        
        %% 4. Transition Matrix for Full System States (MTDC_TM_full)
        % Initialize the full system transition matrix.
        % MTDC_TM_full(i, j) will store the transition rate from system state i to system state j.
        MTDC_TM_full = zeros(num_combinations, num_combinations); 
        
        % Create a map to easily convert state characters ('N', 'P', 'F') to their row/column indices in DCTCP_TM.
        state_indices = containers.Map({'N', 'P', 'F'}, {1, 2, 3});
        
        % Calculate the transition rates between all possible system states.
        % We only consider "single-step" transitions, meaning only one DCTCP changes its state at a time.
        for i = 1:num_combinations
            for j = 1:num_combinations
                if i ~= j % Transitions from a state to itself are handled by diagonal elements later, or implicitly zero for this calculation
                    differences = 0; % Counter for the number of DCTCPs whose states differ between i and j
                    transition_rate = 1; % Initialize transition rate for this (i,j) pair
                    
                    for k = 1:num_MTDC_ends
                        state_i = ternary_table{i, k}; % State of k-th DCTCP in state i
                        state_j = ternary_table{j, k}; % State of k-th DCTCP in state j
                        
                        if ~strcmp(state_i, state_j) % If the k-th DCTCP's state changes
                            differences = differences + 1; % Increment difference count
                            % Get the transition rate for this specific DCTCP's state change
                            idx_i = state_indices(state_i);
                            idx_j = state_indices(state_j);
                            transition_rate = transition_rate * DCTCP_TM(idx_i, idx_j);
                        end
                    end
                    
                    % Only record transitions where exactly one DCTCP changes its state.
                    if differences == 1
                        MTDC_TM_full(i, j) = transition_rate;
                    end
                end
            end
        end
        
        %% 5. Calculate Upward (Lambda Plus) and Downward (Lambda Minus) Transition Rates
        % These rates represent the total transition rate from a given system state
        % to any other state with a higher capacity (lambda_plus) or lower capacity (lambda_minus).
        
        n = size(MTDC_TM_full, 1); % Number of full system states
        lambda_plus_combined = zeros(n, 1);  % Sum of transition rates to higher capacity states
        lambda_minus_combined = zeros(n, 1); % Sum of transition rates to lower capacity states
        
        % Loop through each system state (row 'i') and calculate its lambda_plus and lambda_minus.
        for i = 1:n
            for j = 1:n
                % If state 'j' has a higher capacity than state 'i'
                if capacity_vector(i) < capacity_vector(j)
                    lambda_plus_combined(i) = lambda_plus_combined(i) + MTDC_TM_full(i,j);
                % If state 'j' has a lower capacity than state 'i'
                elseif capacity_vector(i) > capacity_vector(j)
                    lambda_minus_combined(i) = lambda_minus_combined(i) + MTDC_TM_full(i,j);
                end
            end
        end
        
        % Combine all computed information for each full system state into a single table.
        % Column 1: Capacity of the state
        % Column 2: Probability of the state
        % Column 3: Lambda Plus (total rate of moving to a higher capacity state from this state)
        % Column 4: Lambda Minus (total rate of moving to a lower capacity state from this state)
        MTDCCOPT = [capacity_vector, probability_vector, lambda_plus_combined, lambda_minus_combined];
        
        %% 6. Merging Similar Capacity States
        % This section aggregates states with identical capacities into "merged states"
        % to simplify the system representation.
        
        % Combine and find unique normalized capacities from both sending and receiving ends.
        merged_capacities = unique([SDCT_cap_states(:); RDCT_cap_states(:)]);
        merged_capacities = sort(merged_capacities); % Ensure they are sorted in ascending order
        
        % Initialize the merged output table.
        % Column 1: Capacity (actual MW, rounded)
        % Column 2: Merged Probability
        % Column 3: Merged Lambda Plus
        % Column 4: Merged Lambda Minus
        % Column 5: Total departure rate (prob * (lambda_plus + lambda_minus))
        MTDCCOPT_merged = zeros(length(merged_capacities), 5);
        
        % Convert normalized merged capacities to actual MW values and round them.
        MTDCCOPT_merged(:, 1) = round(merged_capacities * MTDC_CAP);
        
        % A cell array to store the original indices (from MTDCCOPT) that map to each merged capacity bin.
        indexes_accumulated = cell(length(merged_capacities), 1); 
        
        % Initialize accumulated frequencies for upward and downward transitions for merged states.
        rounded_freq_plus = zeros(length(merged_capacities), 1);
        rounded_freq_minus = zeros(length(merged_capacities), 1);
        
        % Calculate the merged probabilities for each capacity state.
        % This is the sum of probabilities of all original states that have the same capacity.
        for i = 1:length(merged_capacities)
            MTDCCOPT_merged(i, 2) = sum(probability_vector(capacity_vector == merged_capacities(i)));
        end
        
        % Accumulate raw frequencies (probability * transition rate) into their respective
        % merged capacity bins. Also store which original states fall into each bin.
        for i = 1:size(MTDCCOPT, 1)
            % Find the index of the merged capacity that matches the current original state's capacity.
            [~, bin_idx] = min(abs(merged_capacities - MTDCCOPT(i,1))); % Using min(abs()) for robust float comparison
            
            rounded_freq_plus(bin_idx) = rounded_freq_plus(bin_idx) + MTDCCOPT(i, 2) * MTDCCOPT(i, 3);
            rounded_freq_minus(bin_idx) = rounded_freq_minus(bin_idx) + MTDCCOPT(i, 2) * MTDCCOPT(i, 4);
            
            indexes_accumulated{bin_idx} = [indexes_accumulated{bin_idx}, i]; % Store original index
        end
        
        %% 7. Correction for Transitions Within Merged States
        % When merging states, transitions between original states that end up in the
        % *same* merged capacity state should not contribute to the "net" upward/downward
        % transition rates *between* merged states. This section subtracts these intra-bin transitions.
        
        sum_products_up = zeros(length(merged_capacities), 1);
        sum_products_down = zeros(length(merged_capacities), 1);

        for i = 1:length(merged_capacities)
            indexes = indexes_accumulated{i}; % Get original state indices for the current merged bin 'i'
            
            for j = 1:length(indexes)
                idx_j = indexes(j); % First original state index in the current bin
                for k = 1:length(indexes)
                    idx_k = indexes(k); % Second original state index in the current bin
                    
                    if idx_j ~= idx_k % Exclude self-transitions within the same original state
                        % If transition from idx_j to idx_k is an "upward" transition within the same merged bin
                        if MTDCCOPT(idx_j, 1) < MTDCCOPT(idx_k, 1)
                            sum_products_up(i) = sum_products_up(i) + MTDCCOPT(idx_j, 2) * MTDC_TM_full(idx_j, idx_k);
                        % If transition from idx_j to idx_k is a "downward" transition within the same merged bin
                        elseif MTDCCOPT(idx_j, 1) > MTDCCOPT(idx_k, 1)
                            sum_products_down(i) = sum_products_down(i) + MTDCCOPT(idx_j, 2) * MTDC_TM_full(idx_j, idx_k);
                        end
                    end
                end
            end
        end
        
        % Apply the correction: subtract the intra-bin transition "products" from the accumulated frequencies.
        for i = 1:length(merged_capacities)
            rounded_freq_plus(i) = rounded_freq_plus(i) - sum_products_up(i);
            rounded_freq_minus(i) = rounded_freq_minus(i) - sum_products_down(i) ;
        end
        
        %% 8. Calculate Merged Transition Rates and Total Departure Rates
        
        % Calculate the averaged upward transition rate for each merged state.
        % This is the "net" rate of leaving this merged state to a higher capacity merged state.
        non_zero_prob_idx_plus = MTDCCOPT_merged(:, 2) ~= 0; % Avoid division by zero
        MTDCCOPT_merged(non_zero_prob_idx_plus, 3) = rounded_freq_plus(non_zero_prob_idx_plus) ./ MTDCCOPT_merged(non_zero_prob_idx_plus, 2);
        
        % Calculate the averaged downward transition rate for each merged state.
        % This is the "net" rate of leaving this merged state to a lower capacity merged state.
        non_zero_prob_idx_minus = MTDCCOPT_merged(:, 2) ~= 0; % Avoid division by zero
        MTDCCOPT_merged(non_zero_prob_idx_minus, 4) = rounded_freq_minus(non_zero_prob_idx_minus) ./ MTDCCOPT_merged(non_zero_prob_idx_minus, 2);
        
        % Calculate the total departure rate from each merged state.
        % This is the probability of the state multiplied by its total outward transition rate.
        MTDCCOPT_merged(:, 5) = MTDCCOPT_merged(:, 2) .* (MTDCCOPT_merged(:, 3) + MTDCCOPT_merged(:, 4));
        
        % Flip the merged table upside down. This is typically done to order
        % states from highest capacity to lowest capacity, or vice-versa, depending on convention.
        MTDCCOPT_merged = flipud(MTDCCOPT_merged); 
                
        %% 9. Create the Merged Transition Matrix (MTDC_TM_merged)
        % This matrix represents the transition rates between the merged capacity states.
        
        % Initialize a matrix to store the accumulated frequencies for transitions between merged states.
        freq_rounded = zeros(length(merged_capacities), length(merged_capacities));
        % Initialize the final merged transition matrix.
        MTDC_TM_merged = zeros(length(merged_capacities), length(merged_capacities));
        
        % Calculate the total transition frequency from one merged state to another.
        % This involves summing the products of (original state probability * original transition rate)
        % for all transitions between original states that map to the respective merged states.
        for i = 1:length(merged_capacities) % Loop through source merged states
            indices_i = indexes_accumulated{i}; % Get original state indices belonging to source merged state 'i'
            for j = 1:length(merged_capacities) % Loop through destination merged states
                if i ~= j % Only consider transitions between different merged states
                    indices_j = indexes_accumulated{j}; % Get original state indices belonging to dest merged state 'j'
                    
                    % Sum transitions from all original states in 'indices_i' to all original states in 'indices_j'
                    for idx_i = indices_i
                        for idx_j = indices_j
                            freq_rounded(i, j) = freq_rounded(i, j) + MTDC_TM_full(idx_i, idx_j) * MTDCCOPT(idx_i, 2);
                        end
                    end
                end
            end
        end
        
        % Calculate the final merged transition rates by dividing the accumulated
        % frequencies by the probability of the source merged state.
        % Note: Care must be taken with the indexing of MTDCCOPT_merged due to the `flipud` operation.
        % We need to find the correct row in the flipped `MTDCCOPT_merged` that corresponds to `merged_capacities(i)`.
        for i = 1:length(merged_capacities) % Loop through source merged states
            for j = 1:length(merged_capacities) % Loop through destination merged states
                if i ~= j
                    % Find the row in MTDCCOPT_merged that corresponds to the current source merged capacity.
                    % `round(merged_capacities(i)*MTDC_CAP)` converts the normalized capacity back to its rounded actual value.
                    [~, current_merged_row_idx] = ismember(round(merged_capacities(i)*MTDC_CAP), MTDCCOPT_merged(:,1));
                    
                    if current_merged_row_idx ~= 0 && MTDCCOPT_merged(current_merged_row_idx, 2) ~= 0 % Check if the state exists and has non-zero probability
                     MTDC_TM_merged(i, j) = freq_rounded(i, j) / MTDCCOPT_merged(current_merged_row_idx, 2);
                    else
                     MTDC_TM_merged(i, j) = 0; % If probability is zero, no transitions from this state are possible
                    end
                end
            end
        end

        % After computing MTDC_TM_merged based on the original `merged_capacities` order,
        % it's common practice to flip it to match the order of `MTDCCOPT_merged` if
        % `MTDCCOPT_merged` was `flipud`'d. This ensures consistency between outputs.
        MTDC_TM_merged = flipud(fliplr(MTDC_TM_merged)); % Flip both rows and columns
end

%% Helper Functions (Provided by the user, with added comments and minor adjustments for robustness)

function [capacity_mapping, cap_states] = MTDC_functions_helper(nComponents)
    % MTDC_functions_helper calculates the capacities for all combinations of DCTCP states.
    % It returns a map from state string combinations to their normalized capacities,
    % and a list of all unique normalized capacity values.
    %
    % Inputs:
    %   nComponents - The number of DCTCP components (e.g., nSEDCTCPs or nREDCTCPs).
    %
    % Outputs:
    %   capacity_mapping - A containers.Map object. Keys are string combinations of states
    %                      (e.g., 'NN', 'NPF'), values are their total normalized capacities.
    %   cap_states       - A 1xN array of unique normalized capacity values.

    % Define base normalized capacities for each individual DCTCP state.
    % 'N': 100% (1.0), 'P': 50% (0.5), 'F': 0% (0.0).
    state_capacities_map = containers.Map({'N', 'P', 'F'}, [1, 0.5, 0]);
    
    if nComponents == 0
        % Special case: If there are no components, the combined capacity is 0 (or 1, depending on
        % interpretation in a larger system - here assuming 0 for individual side).
        % An empty string '' represents the state of having no components.
        capacity_mapping = containers.Map();
        capacity_mapping('') = 0; 
        cap_states = 0; % The only capacity state is 0.
        return;
    end

    % Generate all possible combinations of states (e.g., 'NN', 'NP', 'NF', 'PN', ...).
    state_combinations = generate_combinations_helper(nComponents);
    
    % Initialize an array to store the calculated capacities for each combination.
    capacities = zeros(1, length(state_combinations));
    
    % Calculate the total capacity for each combination.
    for i = 1:length(state_combinations)
        combination = state_combinations{i};
        capacities(i) = calculate_combination_capacity_helper(combination, state_capacities_map);
    end
    
    % Create a map where keys are the state string combinations and values are their capacities.
    capacity_mapping = containers.Map(state_combinations, capacities);
    
    % Get unique capacity values and return them.
    cap_states = unique(capacities);
end

function combinations = generate_combinations_helper(n)
    % generate_combinations_helper generates all possible string combinations
    % of DCTCP states for 'n' components.
    %
    % Inputs:
    %   n - The number of components.
    %
    % Outputs:
    %   combinations - A cell array of strings, where each string is a unique
    %                  combination of 'N', 'P', 'F' states.

    states = {'N', 'P', 'F'}; % Possible states for a single DCTCP
    num_base_states = length(states); % Number of possible states (3)
    
    if n == 0
        combinations = {''}; % If no components, the only "combination" is an empty string
        return;
    end

    % Pre-allocate cell array for efficiency. Total combinations = (num_states)^n.
    combinations = cell(1, num_base_states^n);
    idx = 1; % Index for placing combinations into the cell array
    
    % Loop from 0 to (3^n - 1) to generate all combinations.
    % This uses a base-3 counting method.
    for i = 0:(num_base_states^n - 1)
        combination = ''; % Initialize an empty string for the current combination
        num = i;          % Temporary variable for the current count
        
        % Build the combination string by extracting digits in base 3.
        % The states are appended from right to left to match the original code's logic.
        for j = 1:n
            state_idx = mod(num, num_base_states) + 1; % Get the remainder (0, 1, or 2), add 1 for 1-based indexing
            combination = strcat(states{state_idx}, combination); % Prepend the state character
            num = floor(num / num_base_states); % Move to the next "digit"
        end
        combinations{idx} = combination; % Store the generated combination
        idx = idx + 1;
    end
end

function capacity = calculate_combination_capacity_helper(combination, state_capacities_map)
    % calculate_combination_capacity_helper calculates the average capacity
    % for a given string combination of DCTCP states.
    %
    % Inputs:
    %   combination        - A string representing the combined states (e.g., 'NPF').
    %   state_capacities_map - A containers.Map mapping single state chars ('N','P','F') to their capacities.
    %
    % Outputs:
    %   capacity - The average normalized capacity for the given combination.

    capacity = 0; % Initialize total capacity for this combination
    
    if isempty(combination)
        % For an empty combination (e.g., when n=0), capacity is 0.
        return;
    end
    
    % Iterate through each character in the combination string.
    for i = 1:length(combination)
        state_char = combination(i); % Get the character representing a single DCTCP's state
        capacity = capacity + state_capacities_map(state_char); % Add its capacity
    end
    
    % The combined capacity is the average of individual DCTCP capacities.
    capacity = capacity / length(combination);
end