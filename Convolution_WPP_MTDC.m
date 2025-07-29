function [CONVOLVED_SCOPT_merged, CONVOLVED_SCOPT_TM_merged] = Convolution_WPP_MTDC(FCOPT_TM_rounded, MTDC_TM_merged, rounded_FCOPT, MTDCCOPT_merged)
% CONVOLUTION_WPP_MTDC performs the convolution of WPP and MTDC probabilistic models.
% This function combines the state probabilities and transition matrices of
% a Wind Power Plant (WPP) model and an MTDC model to create a convolved
% system model. It then merges similar states in the convolved model.
%
% Inputs:
%   FCOPT_TM_rounded  - Rounded transition matrix for the WPP states.
%   MTDC_TM_merged    - Merged transition matrix for the MTDC system states.
%   rounded_FCOPT     - Table of WPP states (capacity, probability).
%                        Column 1: Capacity, Column 2: Probability.
%   MTDCCOPT_merged   - Merged table of MTDC states (capacity, probability, rates).
%                        Column 1: Capacity, Column 2: Probability,
%                        Column 3: Upward rate, Column 4: Downward rate, Column 5: Total departure rate.
%
% Outputs:
%   CONVOLVED_SCOPT_merged    - Merged convolved system table.
%                               Column 1: Capacity, Column 2: Probability,
%                               Column 3: Upward rate, Column 4: Downward rate, Column 5: Total departure rate.
%   CONVOLVED_SCOPT_TM_merged - Merged convolved system transition matrix.

        % Convolved WPP+MTDC transition matrix
        
        % blockSize corresponds to the number of states in the WPP model (rows/cols of FCOPT_TM_rounded)
        blockSize = size(FCOPT_TM_rounded, 1);
        % totalSize is the dimension of the convolved system's transition matrix.
        % It's the product of the number of WPP states and the number of MTDC states.
        totalSize = blockSize * size(MTDC_TM_merged, 1);
        
        % Initialize the large convolved transition matrix with zeros.
        CONVOLVED_SCOPT_TM = zeros(totalSize, totalSize);
        
        % Loop to construct the CONVOLVED_SCOPT_TM by placing blocks.
        % Each block represents transitions within WPP states, conditioned on an MTDC state.
        % Transitions between blocks represent MTDC state changes.
        for i = 1:size(MTDC_TM_merged, 1) % Iterate through MTDC states (rows of MTDC_TM_merged)
            % Calculate the starting row and column for the current block in the large matrix.
            % Each MTDC state corresponds to a 'blockSize' x 'blockSize' block.
            rowStart = (i-1) * blockSize + 1;
            colStart = rowStart;
            
            % Place the FCOPT_TM_rounded block (WPP internal transitions) along the diagonal.
            % This represents transitions within the WPP states while the MTDC state 'i' remains constant.
            CONVOLVED_SCOPT_TM(rowStart:rowStart+blockSize-1, colStart:colStart+blockSize-1) = FCOPT_TM_rounded;
            
            % Assign transition rates corresponding to MTDC state changes to higher states (i+1, i+2).
            % These are placed off-diagonal, indicating transitions from MTDC state 'i' to 'i+1' or 'i+2',
            % while maintaining the relative WPP state within the block.
            
            % Transition to MTDC state i+1
            if i < size(MTDC_TM_merged, 1) % Ensure there's a next MTDC state
                % Loop through each internal WPP state within the current block
                for j = 1:blockSize
                    % Check if the target column index is within the bounds of CONVOLVED_SCOPT_TM
                    if colStart + blockSize + (j-1) <= totalSize 
                        % Assign the transition rate from MTDC state i to i+1.
                        % This rate applies to all WPP state transitions within the (i, i+1) block.
                        CONVOLVED_SCOPT_TM(rowStart+j-1, colStart+blockSize+j-1) = MTDC_TM_merged(i, i+1);
                    end
                    % Transition to MTDC state i+2
                    if i < size(MTDC_TM_merged, 1) - 1 % Ensure there's an MTDC state i+2
                         if colStart + 2*blockSize + (j-1) <= totalSize
                            CONVOLVED_SCOPT_TM(rowStart+j-1, colStart+2*blockSize+j-1) = MTDC_TM_merged(i, i+2);
                        end
                    end
                end
            end
            
            % Assign transition rates corresponding to MTDC state changes to lower states (i-1, i-2).
            % These are also placed off-diagonal, indicating transitions from MTDC state 'i' to 'i-1' or 'i-2'.
            
            % Transition to MTDC state i-1
            if i > 1 % Ensure there's a previous MTDC state
                % Loop through each internal WPP state within the current block
                for j = 1:blockSize
                    % Check if the target column index is within the bounds of CONVOLVED_SCOPT_TM
                    if colStart - blockSize + (j-1) > 0 % Ensure colStart-blockSize+j-1 is positive and valid
                        CONVOLVED_SCOPT_TM(rowStart+j-1, colStart-blockSize+j-1) = MTDC_TM_merged(i, i-1);
                    end
                    % Transition to MTDC state i-2
                    if i > 2 % Ensure there's an MTDC state i-2
                        if colStart - 2*blockSize + (j-1) > 0
                            CONVOLVED_SCOPT_TM(rowStart+j-1, colStart-2*blockSize+j-1) = MTDC_TM_merged(i, i-2);
                        end
                    end
                end
            end
        end
        
        %% 2. Calculate Combined Probabilities (Convolution of Capacities and Probabilities)
        
        % Initialize empty arrays to store the combined capacities and probabilities.
        P_combined = [];       % Stores the minimum capacity (convolved capacity)
        prob_combined = [];    % Stores the combined probability
        
        % Iterate through all combinations of MTDC states and WPP states.
        % The combined system capacity is the minimum of the WPP capacity and the MTDC capacity.
        % The combined probability is the product of their individual probabilities (assuming independence).
        for i = 1:size(MTDCCOPT_merged, 1) % Loop through MTDC merged states
            for j = 1:size(rounded_FCOPT, 1) % Loop through WPP rounded states
                % The convolved output capacity is the minimum of the WPP capacity and MTDC capacity.
                P_new = min(rounded_FCOPT(j, 1), MTDCCOPT_merged(i, 1));
                % The combined probability is the product of their individual probabilities.
                prob_new = MTDCCOPT_merged(i, 2) * rounded_FCOPT(j, 2);
        
                % Append the new combined capacity and probability to the lists.
                P_combined = [P_combined; P_new];
                prob_combined = [prob_combined; prob_new];
            end
        end
        
        %% 3. Calculate Upward (Lambda Plus) and Downward (Lambda Minus) Transition Rates for Convolved States
        % These rates represent the total transition rate from a given convolved state
        % to any other convolved state with a higher capacity (lambda_plus) or lower capacity (lambda_minus).
        
        n = size(CONVOLVED_SCOPT_TM, 1); % Number of convolved system states
        lambda_plus_combined = zeros(n, 1);  % Sum of transition rates to higher capacity states
        lambda_minus_combined = zeros(n, 1); % Sum of transition rates to lower capacity states
        
        % Loop through each convolved system state (row 'i') and calculate its lambda_plus and lambda_minus.
        for i = 1:n
            for j = 1:n
                % If state 'j' has a higher capacity than state 'i'
                if P_combined(i) < P_combined(j)
                    lambda_plus_combined(i) = lambda_plus_combined(i) + CONVOLVED_SCOPT_TM(i,j);
                % If state 'j' has a lower capacity than state 'i'
                elseif P_combined(i) > P_combined(j)
                    lambda_minus_combined(i) = lambda_minus_combined(i) + CONVOLVED_SCOPT_TM(i,j);
                end
            end
        end
        
        % Create the final combined system table (before merging identical capacities).
        % Column 1: Convolved Capacity
        % Column 2: Convolved Probability
        % Column 3: Lambda Plus (total rate of moving to a higher capacity convolved state)
        % Column 4: Lambda Minus (total rate of moving to a lower capacity convolved state)
        CONVOLVED_SCOPT = [P_combined, prob_combined, lambda_plus_combined, lambda_minus_combined];
        
        %% 4. Merging Identical Capacity States in the Convolved Table
        % This section aggregates convolved states with identical capacities into "merged states"
        % to simplify the system representation. This is similar to the merging process in MTDC_Probabilistic_Modelling.
        
        % Define the unique capacity bins from the convolved system.
        bins = unique(CONVOLVED_SCOPT(:,1));
        bins = sort(bins); % Ensure sorted order for consistency
        
        % Initialize the merged convolved table.
        % Column 1: Merged Capacity
        % Column 2: Merged Probability
        % Column 3: Merged Upward rate
        % Column 4: Merged Downward rate
        % Column 5: Total departure rate
        CONVOLVED_SCOPT_merged = zeros(length(bins), 5);
        CONVOLVED_SCOPT_merged(:, 1) = bins;
        
        % A cell array to store the original convolved state indices that map to each merged capacity bin.
        indexes_accumulated = cell(length(bins), 1); 
        
        % Initialize accumulated frequencies for upward and downward transitions for merged states.
        rounded_freq_plus = zeros(length(bins), 1);
        rounded_freq_minus = zeros(length(bins), 1);
        
        % Calculate the merged probabilities for each capacity bin.
        % This is the sum of probabilities of all original convolved states that have the same capacity.
        for i = 1:length(bins)
            CONVOLVED_SCOPT_merged(i, 2) = sum(CONVOLVED_SCOPT(CONVOLVED_SCOPT(:,1) == bins(i), 2));
        end
        
        % Accumulate raw frequencies (probability * transition rate) into their respective
        % merged capacity bins. Also store which original convolved states fall into each bin.
        for i = 1:size(CONVOLVED_SCOPT, 1)
            % Find the index of the merged capacity bin that matches the current original convolved state's capacity.
            bin_idx = find(bins == CONVOLVED_SCOPT(i,1)); % Assumes exact match for binning
            
            rounded_freq_plus(bin_idx) = rounded_freq_plus(bin_idx) + CONVOLVED_SCOPT(i, 2) * CONVOLVED_SCOPT(i, 3);
            rounded_freq_minus(bin_idx) = rounded_freq_minus(bin_idx) + CONVOLVED_SCOPT(i, 2) * CONVOLVED_SCOPT(i, 4);
            
            indexes_accumulated{bin_idx} = [indexes_accumulated{bin_idx}, i]; % Store original index
        end
        
        %% 5. Correction for Transitions Within Merged Convolved States
        % Similar to the MTDC function, this corrects for transitions between original
        % convolved states that fall into the same merged capacity bin, ensuring they
        % don't contribute to net inter-bin transition rates.
        
        sum_products_up = zeros(length(bins), 1);
        sum_products_down = zeros(length(bins), 1);
        
        for i = 1:length(bins)
            indexes = indexes_accumulated{i}; % Get original convolved state indices for the current merged bin 'i'
            
            for j = 1:length(indexes)
                idx_j = indexes(j); % First original convolved state index in the current bin
                for k = 1:length(indexes)
                    idx_k = indexes(k); % Second original convolved state index in the current bin
                    
                    if idx_j ~= idx_k % Exclude self-transitions within the same original state
                        % If transition from idx_j to idx_k is an "upward" transition within the same merged bin
                        if CONVOLVED_SCOPT(idx_j, 1) < CONVOLVED_SCOPT(idx_k, 1)
                            sum_products_up(i) = sum_products_up(i) + CONVOLVED_SCOPT(idx_j, 2) * CONVOLVED_SCOPT_TM(idx_j, idx_k);
                        % If transition from idx_j to idx_k is a "downward" transition within the same merged bin
                        elseif CONVOLVED_SCOPT(idx_j, 1) > CONVOLVED_SCOPT(idx_k, 1)
                            sum_products_down(i) = sum_products_down(i) + CONVOLVED_SCOPT(idx_j, 2) * CONVOLVED_SCOPT_TM(idx_j, idx_k);
                        end
                    end
                end
            end
        end
        
        % Apply the correction to the accumulated frequencies.
        for i = 1:length(bins) % Loop through merged bins, not MTDCCOPT_merged
            rounded_freq_plus(i) = rounded_freq_plus(i) - sum_products_up(i);
            rounded_freq_minus(i) = rounded_freq_minus(i) - sum_products_down(i);
        end
        
        %% 6. Calculate Merged Transition Rates and Total Departure Rates for Convolved System
        
        % Calculate the averaged upward transition rate for each merged convolved state.
        non_zero_prob_idx_plus = CONVOLVED_SCOPT_merged(:, 2) ~= 0; % Avoid division by zero
        CONVOLVED_SCOPT_merged(non_zero_prob_idx_plus, 3) = rounded_freq_plus(non_zero_prob_idx_plus) ./ CONVOLVED_SCOPT_merged(non_zero_prob_idx_plus, 2);
        
        % Calculate the averaged downward transition rate for each merged convolved state.
        non_zero_prob_idx_minus = CONVOLVED_SCOPT_merged(:, 2) ~= 0; % Avoid division by zero
        CONVOLVED_SCOPT_merged(non_zero_prob_idx_minus, 4) = rounded_freq_minus(non_zero_prob_idx_minus) ./ CONVOLVED_SCOPT_merged(non_zero_prob_idx_minus, 2);
        
        % Calculate the total departure rate from each merged convolved state.
        CONVOLVED_SCOPT_merged(:, 5) = CONVOLVED_SCOPT_merged(:, 2) .* (CONVOLVED_SCOPT_merged(:, 3) + CONVOLVED_SCOPT_merged(:, 4));
        
        % Flip the merged table upside down.
        CONVOLVED_SCOPT_merged = flipud(CONVOLVED_SCOPT_merged);
        
        %% 7. Create the Merged Convolved Transition Matrix (CONVOLVED_SCOPT_TM_merged)
        % This matrix represents the transition rates between the merged convolved states.
        
        % Initialize matrices for accumulated frequencies and the final merged transition matrix.
        freq_rounded = zeros(length(bins), length(bins));
        CONVOLVED_SCOPT_TM_merged = zeros(length(bins), length(bins));
        
        % Calculate the total transition frequency from one merged convolved state to another.
        for i = 1:length(bins) % Loop through source merged states
            indices_i = indexes_accumulated{i}; % Original convolved state indices for source merged state 'i'
            for j = 1:length(bins) % Loop through destination merged states
                if i ~= j % Only consider transitions between different merged states
                    indices_j = indexes_accumulated{j}; % Original convolved state indices for dest merged state 'j'
                    
                    % Sum transitions from all original convolved states in 'indices_i' to
                    % all original convolved states in 'indices_j'.
                    for idx_i = indices_i
                        for idx_j = indices_j
                            freq_rounded(i, j) = freq_rounded(i, j) + CONVOLVED_SCOPT_TM(idx_i, idx_j) * CONVOLVED_SCOPT(idx_i, 2);
                        end
                    end
                end
            end
        end
        
        % Calculate the final merged convolved transition rates.
        % Divide the accumulated frequencies by the probability of the source merged state.
        % Note: Similar to the previous function, handling `flipud` for `CONVOLVED_SCOPT_merged` is crucial here.
        % We need to find the correct row in the flipped `CONVOLVED_SCOPT_merged` corresponding to `bins(i)`.
        for i = 1:length(bins) % Loop through source merged states
            for j = 1:length(bins) % Loop through destination merged states
                if i ~= j
                    % Find the row in CONVOLVED_SCOPT_merged that corresponds to the current source merged capacity.
                    % Using `ismember` or `find` to handle cases where `flipud` reordered the table.
                    [~, current_merged_row_idx] = ismember(bins(i), CONVOLVED_SCOPT_merged(:,1));
                    
                    if current_merged_row_idx ~= 0 && CONVOLVED_SCOPT_merged(current_merged_row_idx, 2) ~= 0 % Check for non-zero probability
                        CONVOLVED_SCOPT_TM_merged(i, j) = freq_rounded(i, j) / CONVOLVED_SCOPT_merged(current_merged_row_idx, 2);
                    else
                        CONVOLVED_SCOPT_TM_merged(i, j) = 0; % If probability is zero, transition rate is zero
                    end
                end
            end
        end

        % Flip the merged transition matrix to match the order of `CONVOLVED_SCOPT_merged`.
        CONVOLVED_SCOPT_TM_merged = flipud(fliplr(CONVOLVED_SCOPT_TM_merged));
end