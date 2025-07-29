function [rounded_FCOPT, FCOPT_TM_rounded, WT_States] = Wind_Farm_Probabilistic_Modelling(turbine, WPP_CAP, offshore_distance, num_wind_states, num_rounding_states)
%Wind_Farm_Probabilistic_Modelling Calculates reliability metrics for an offshore wind power plant.
%
%   [rounded_FCOPT, FCOPT_TM_rounded] = Wind_Farm_Probabilistic_Modelling(turbine, WPP_CAP, offshore_distance, num_wind_states, num_rounding_states);
%   This function determines the reliability characteristics of a wind power plant
%   by combining turbine outage probabilities with wind speed variability. It
%   outputs two key matrices: a 'Farm Capacity Output Probability Table'
%   (rounded_FCOPT) and a corresponding 'Rounded System Outage Probability
%   Transition Matrix' (FCOPT_TM_rounded).
%
%   Required Inputs:
%   - turbine (char array): Specifies the turbine model. Accepted values are
%     'IEA' or 'DTU'. This input dictates the turbine's rated power,
%     failure rate, and wind-power curve characteristics.
%   - WPP_CAP (double): The total installed capacity of the wind power plant in MW.
%     This is used to determine the total number of turbines required.
%   - offshore_distance (double): The distance of the offshore wind farm from shore in km.
%     This distance influences the Forced Outage Rate (FOR) of the turbines.
%   - num_wind_states (double): The number of discrete wind power output states
%     to model for the wind speed probability distribution. Higher values offer
%     finer granularity in wind power modeling.
%   - num_rounding_states (double or char array): Controls the number of
%     capacity states in the final 'Farm Capacity Output Probability Table'.
%     If a numeric value, it sets the number of equally spaced output capacity
%     bins. If 'exact', all unique combined capacity values are preserved.
%
%   Outputs:
%   - rounded_FCOPT (matrix): A table representing the aggregated
%     Farm Capacity Output Probability. Its columns typically include:
%     [Output Capacity (MW), Probability, Upward Transition Rate,
%     Downward Transition Rate, Frequency]. It is ordered from highest
%     to lowest output capacity.
%   - FCOPT_TM_rounded (matrix): A square matrix detailing the transition rates
%     between the rounded farm capacity output states.

%% Turbine and Plant Initialization
% This section defines the specific parameters for the chosen turbine type
% and calculates the total number of turbines required for the wind power plant.

        if strcmp(turbine,'IEA')
            wT= 15; % IEA Turbine rated power
            lambda_wT=38.5396557297363/365; % IEA Equivalent failure rate in [occ/day] Reference: James Carroll et al. DOI:10.1002/we.2011
        elseif strcmp(turbine,'DTU')
            wT=10; % DTU Turbine rated power
            lambda_wT=39.86052784/365; % DTU Equivalent failure rate in [occ/day] Reference: James Carroll et al. DOI:10.1002/we.2011
        else
            error('Invalid turbine type. Please choose ''IEA'' or ''DTU''.');
        end
        
        nT=round(WPP_CAP/wT); % Number of turbines required
        units = wT*ones(nT,1);  % Capacities in MW
        
        failure_rate = lambda_wT*ones(nT,1); % Units failure rates in [occ/day]
        
%% Forced Outage Rate (FOR) Correction based on Offshore Distance
% This section corrects the turbine's Forced Outage Rate (FOR) to account
% for the offshore distance, as maintenance and repair times can increase
% with distance from shore. Different formulas are applied based on distance ranges.

        FOR_corrected = 0; 

        if strcmp(turbine,'IEA')    
            if offshore_distance <= 40 && offshore_distance >= 10    
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+(offshore_distance-10)*0.004/30; % Corrected FOR probabilities to distance 10-40 km
            elseif offshore_distance <= 50 && offshore_distance > 40
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+30*0.004/30+(offshore_distance-40)*0.0005; % Corrected FOR probabilities to distance 40-50 km
            elseif offshore_distance <= 70 && offshore_distance > 50
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+30*0.004/30+10*0.0005+(offshore_distance-50)*0.00015; % Corrected FOR probabilities to distance 50-70 km
            elseif offshore_distance <= 80 && offshore_distance > 70
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+30*0.004/30+10*0.0005+20*0.00015+(offshore_distance-70)*0.0012; % Corrected FOR probabilities to distance 70-80 km
            elseif offshore_distance <= 90 && offshore_distance > 80
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+30*0.004/30+10*0.0005+20*0.00015+10*0.0012+(offshore_distance-80)*0.0078; % Corrected FOR probabilities to distance 80-90 km
            elseif offshore_distance > 90
                FOR_corrected = (38.5396557297363/(38.5396557297363+817.190322341658))+30*0.004/30+10*0.0005+20*0.00015+10*0.0012+10*0.0078+(offshore_distance-90)*0.0008; % Corrected FOR probabilities to distance >90 km
            else % Added error handling for distances outside the specified ranges
                error('Offshore distance for IEA turbine not in specified range (>=10 km).');
            end
            
        elseif strcmp(turbine,'DTU')
            if offshore_distance <= 40 && offshore_distance >= 10    
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+(offshore_distance-10)*0.00015; % Corrected FOR probabilities to distance 10-40 km
            elseif offshore_distance <= 50 && offshore_distance > 40
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+30*0.00015+(offshore_distance-40)*0.0007; % Corrected FOR probabilities to distance 40-50 km
            elseif offshore_distance <= 70 && offshore_distance > 50
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+30*0.00015+10*0.0007+(offshore_distance-50)*0.00025; % Corrected FOR probabilities to distance 50-70 km
            elseif offshore_distance <= 80 && offshore_distance > 70
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+30*0.00015+10*0.0007+20*0.00025+(offshore_distance-70)*0.0013; % Corrected FOR probabilities to distance 70-80 km
            elseif offshore_distance <= 90 && offshore_distance > 80
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+30*0.00015+10*0.0007+20*0.00025+10*0.0013+(offshore_distance-80)*0.0086;% Corrected FOR probabilities to distance 80-90 km
            elseif offshore_distance > 90
                FOR_corrected = (39.86052784/(39.86052784+601.550856174534))+30*0.00015+10*0.0007+20*0.00025+10*0.0013+10*0.0086+(offshore_distance-90)*0.001;% Corrected FOR probabilities to distance >90 km
            else % Added error handling for distances outside the specified ranges
                error('Offshore distance for DTU turbine not in specified range (>=10 km).');
            end
        else
            error('Invalid turbine type. Please choose ''IEA'' or ''DTU''.'); % Redundant with initial check, but harmless.
        end
        
        probabilities = FOR_corrected*ones(nT,1);  % FOR probabilities
        
        mu_wT = (lambda_wT*(1-FOR_corrected)/FOR_corrected); % Corrected Equivalent repair rate in [occ/day]
        repair_rate = mu_wT*ones(nT,1); % Units repair rate in [occ/day]
%% Farm Capacity Outage Probability Table (FCOPT) Construction
% This segment builds the Farm Capacity Outage Probability Table (FCOPT) by iteratively adding each turbine and calculating the probabilities of various outage states, along with their associated upward ($\lambda^+$) and downward ($\lambda^-$) transition rates. This table initially represents the **outage capacity** of the wind farm due to turbine failures alone, from 0 MW out (all turbines available) to max capacity out (all turbines failed).
        
        %Initialization for array sizes
        max_capacity = sum(units); % Define the maximum capacity to initialize arrays
        com_cap = 0; %Initialize an accumulating capacity variable for each i unit added
        
        % Initialize the COPT array and transition rates
        FCOPT = []; % Final System Outage probability table from 0 out to max_capacity (all units out) and all its possible intermediate combinations
        TCOPT = zeros(max_capacity + 1, 1);% Recursive outage probability table from 0 out to max_capacity 
        flp_COPT = zeros(max_capacity + 1, 1); % Flipped recursive outage probability table from 0 out to max_capacity
        lambda_plus = zeros(max_capacity + 1, 1);  % Upward transitions
        lambda_minus = zeros(max_capacity + 1, 1);  % Downward transitions
        
        % Start with the system fully available
        TCOPT(1) = 1;
        
        % Loop through each unit
        for i = 1:length(units)
            % Temporary copy to update probabilities and transition rates
            temp_COPT = zeros(max_capacity + 1, 1);
            temp_COPTU = zeros(max_capacity + 1, 1);
            temp_COPTD = zeros(max_capacity + 1, 1);
            temp_lambda_plus = zeros(max_capacity + 1, 1);
            temp_lambda_minus = zeros(max_capacity + 1, 1);
            lmdp = zeros(max_capacity + 1, 1);
            lmdm = zeros(max_capacity + 1, 1);
        
            lmdm(1:com_cap+1)=flipud(lambda_minus(1:com_cap+1));% Flipped recursive decreasing transition rates
        
            for j = 0:max_capacity
                if TCOPT(j + 1) > 0
                    % Probability unit is up
                    if j + units(i) <= max_capacity
                        temp_COPT(j + units(i) + 1) = temp_COPT(j + units(i) + 1) + TCOPT(j + 1) * (1 - probabilities(i));
                        
                    end
                    % Probability unit is down
                    temp_COPT(j + 1) = temp_COPT(j + 1) + TCOPT(j + 1) * probabilities(i);
                    
                    if i~=1
                        temp_COPTU(j + 1) = TCOPT(j + 1) * probabilities(i); % Calculate p'(X-C)*U (Col 3)
                        temp_lambda_plus(j + 1) = temp_COPTU(j + 1)*(lambda_plus(j + 1)+repair_rate(i)); % Calculate p'(X-C)*U*(λ+(X-C)+μ(i)) (Col 6)
                    
                        temp_COPTD(com_cap + units(i) - j + 1) = flp_COPT(j + 1) * (1 - probabilities(i)); % Calculate p'(X)*(1-U) (Col 3)
                        temp_lambda_minus(com_cap + units(i) - j + 1) = temp_COPTD(com_cap + units(i) - j + 1)*(lmdm(j + 1)+ failure_rate(i)); % Calculate p'(X)*(1-U)*(λ-(X)+λ(i)) (Col 9)
                    end
                end
            end
        
            % Update COPT, flipped COPT, and transition rates with the temporary values
            TCOPT = temp_COPT;
            flp_COPT(1:com_cap+units(i)+1)=flipud(temp_COPT(1:com_cap+units(i)+1));
        
            if i==1
                lambda_plus(1) = repair_rate(i); % Initialization of first unit repair rate
                lambda_minus(units(i) + 1) = failure_rate(i); % Initialization of first unit failure rate
            else
        
                lmdp(units(i)+1:com_cap+units(i)+1)=lambda_plus(1:com_cap+1); % Reordered recursive increasing transition rates
                
                temp_lambda_plus = temp_lambda_plus + temp_COPTD.*lmdp; % Calculate p'(X-C)*U*(λ+(X-C)+μ(i)) + p'(X)*(1-U)*λ+(X) (Col 7)
                temp_lambda_minus = temp_lambda_minus + temp_COPTU.*lambda_minus; % Calculate p'(X)*(1-U)*(λ-(X)+λ(i)) +p'(X-C)*U*λ-(X-C) (Col 11)
                lambda_plus = temp_lambda_plus./TCOPT; lambda_plus(isnan(lambda_plus)) = 0;  % Update new λ+(X) and replace all NaN (where there is division by 0)
                lambda_minus = temp_lambda_minus./TCOPT; lambda_minus(isnan(lambda_minus)) = 0; % Update new λ-(X) and replace all NaN (where there is division by 0)
            end
            com_cap = com_cap + units(i); % Accumulate the capacity of the added unit
        end
        % Create the vector X with values from 0 to max_capacity in steps of wT
        X = 0:wT:max_capacity;
        
        % Calculate the required columns for FCOPT
        TCOPT_values = TCOPT(X + 1);
        lambda_plus_values = lambda_plus(X + 1);
        lambda_minus_values = lambda_minus(X + 1);
        frequency_values = TCOPT_values .* (lambda_minus_values + lambda_plus_values);
        
        % Construct the FCOPT matrix
        FCOPT = [X', TCOPT_values, lambda_plus_values, lambda_minus_values, frequency_values];
        
        FCOPT= flipud(FCOPT); % Flip values to start from 0 out to max_capacity (all units out) and all its possible intermediate combinations
%% Wind Speed Probability and Transition Modeling
% This section processes historical wind speed data to create a probability distribution of wind power output states and their corresponding transition rates. This models the variability introduced by wind conditions.
        
        options = [20, 30, 50, 60, 70, 80, 90, 100]; 
        
        % Rounding according to distance
        [~, idx] = min(abs(options - offshore_distance)); 
        rounded_distance = options(idx); 
        
        % Load windspeed data files into a struct
        dataFiles = {'KM_20.mat', 'KM_30.mat', 'KM_50.mat', 'KM_60.mat', 'KM_70.mat', 'KM_80.mat', 'KM_90.mat', 'KM_100.mat'};
        dataStruct = struct();
        
        for i = 1:length(dataFiles)
            % Ensure the .mat files exist in the MATLAB path
            if exist(dataFiles{i}, 'file') == 2
                data = load(dataFiles{i});
                fieldName = sprintf('KM_%d', options(i));
                % Assuming each .mat file contains a single variable, which is the matrix
                temp_field = fieldnames(data);
                dataStruct.(fieldName) = data.(temp_field{1});
            else
                error('Missing data file: %s. Please ensure all KM_XX.mat files are in the MATLAB path.', dataFiles{i});
            end
        end
        
        % Select the appropriate wind speed column based on turbine type and distance
        fieldName = sprintf('KM_%d', rounded_distance);
        if isfield(dataStruct, fieldName)
            windData = dataStruct.(fieldName);
        else
            error('Wind data for rounded_distance %d not found in dataStruct. Check dataFiles and offshore_distance input.', rounded_distance);
        end

        % Select the appropriate wind speed column based on turbine type and distance
        rho=1.225; % Selected air density
        Cp = 0; Area_wT = 0; % Initialize to avoid potential unassigned warnings
        if strcmp(turbine, 'IEA')
            w_speed = windData(:, 2);
            w_speed_cutin = 3; w_speed_cutout=25; Cp=0.455811125095624; Area_wT=pi()*120^2;
        elseif strcmp(turbine, 'DTU')
            w_speed = windData(:, 1);
            w_speed_cutin = 4; w_speed_cutout=25; Cp=0.441353713697127; Area_wT=pi()*89.15^2;
        end

        % Generate rounded states vector
        rounded_p_states = wT*10^6*linspace(0, 1, num_wind_states);
        rounded_w_states = (rounded_p_states*2/(rho*Cp*Area_wT)).^(1/3);

        % Define nested function for roundToNearest to be accessible within this function
        function val_rounded = roundToNearest(value, array_to_round_to)
            [~, idx_nearest] = min(abs(array_to_round_to - value));
            val_rounded = array_to_round_to(idx_nearest);
        end
        
        w_speed_rounded = zeros(size(w_speed)); % Pre-allocate for efficiency
        if strcmp(turbine, 'IEA')
            for i = 1:numel(w_speed)
                if w_speed(i) < w_speed_cutin || w_speed(i) > w_speed_cutout
                    w_speed_rounded(i) = 0;
                else
                    w_speed_rounded(i) = roundToNearest(w_speed(i), rounded_w_states(2:end));
                end
            end
        elseif strcmp(turbine, 'DTU')    
            for i = 1:numel(w_speed)
                if w_speed(i) < w_speed_cutin || w_speed(i) > w_speed_cutout
                    w_speed_rounded(i) = 0;
                else
                    w_speed_rounded(i) = roundToNearest(w_speed(i), rounded_w_states(2:end));
                end
            end
        end
        
        % Count the occurrences of each rounded state
        occurrences = zeros(size(rounded_w_states));
        total_points = numel(w_speed_rounded);
        
        for i = 1:numel(rounded_w_states)
            occurrences(i) = sum(w_speed_rounded == rounded_w_states(i));
        end
        
        % Calculate the probabilities
        w_probabilities = occurrences / total_points;
        w_probabilities = flip(w_probabilities);
        
        % Initialize transition rate matrix
        w_transition_matrix = zeros(num_wind_states);
        
        % Count transitions from state i to state k and divide by occurrences of state i
        for i = 1:numel(rounded_w_states)
            state_i = rounded_w_states(i);
            occurrences_i = occurrences(i);
            
            % Find transitions from state i to state k
            for k = 1:numel(rounded_w_states)
                state_k = rounded_w_states(k);
                
                % Count transitions from state i to state k
                transitions_ik = sum(w_speed_rounded(1:end-1) == state_i & w_speed_rounded(2:end) == state_k);
                
                % Calculate transition rate from state i to state k
                if occurrences_i > 0 && state_i ~= state_k % Exclude staying in the same state
                    transition_rate = transitions_ik / occurrences_i;
                else
                    transition_rate = 0; % Avoid division by zero and set diagonal elements to 0
                end
                
                % Update transition matrix
                w_transition_matrix(i, k) = transition_rate;
            end
        end
        
        w_transition_matrix=w_transition_matrix*24; % Multiply by hours in a day to get transitions in [occ\day]
        
        % Mirror the matrix both horizontally and vertically
        w_transition_matrix = flip(flip(w_transition_matrix, 2),1);
        WT_States = zeros(num_wind_states,5); % Initialize with correct dimensions
        WT_States(:,1) = flip(rounded_p_states')/1000000; WT_States(:,2) = w_probabilities';
        % Loop to build λ+ and λ- dynamically from transition matrix
        for i = 1:num_wind_states
            if i > 1
                WT_States(i,3) = sum(w_transition_matrix(i, 1:i-1)); 
            else
                WT_States(i,3) = 0;  % First row: no upward transition
            end
            if i < num_wind_states
                WT_States(i,4) = sum(w_transition_matrix(i, i+1:end)); 
            else
                WT_States(i,4) = 0;  % Last row: no downward transition
            end
        end
        WT_States(:,5) =  WT_States(:,2).*(WT_States(:,3) + WT_States(:,4));
%% Combined Wind Farm and Wind-Dependent Probability Table (FCOPT)
% This section integrates the turbine failure probabilities with the wind speed probabilities to create a comprehensive system outage/output probability table. It also constructs a large transition matrix (`FCOPT_TM`) that captures all possible transitions between combined turbine outage and wind power states.
        
        % WPP Failure Transition matrix
        % Create an nT x nT matrix initialized to zero
        WPP_transition_matrix = zeros(nT+1, nT+1);
        
        % Set values one position to the right of the diagonal
        for i = 1:nT
            WPP_transition_matrix(i, i+1) = (nT - i + 1) * lambda_wT; 
        end
        
        % Set values one position to the left of the diagonal
        for i = 2:nT+1
            WPP_transition_matrix(i, i-1) = (i-1) * mu_wT;
        end
        
        % Size of the WT_States matrix to determine the block size
        blockSize = size(WT_States, 1);
        totalSize = blockSize * (nT + 1);
        
        % Initialize the large transition matrix
        FCOPT_TM = zeros(totalSize, totalSize);
        
        % Loop to place blocks in the large matrix
        for i = 1:nT+1
            rowStart = (i-1) * blockSize + 1;
            colStart = rowStart;
            
            % Place the transition_matrix block
            FCOPT_TM(rowStart:rowStart+blockSize-1, colStart:colStart+blockSize-1) = w_transition_matrix;
            
            % Assign the right side diagonal values if not in the last block
            if i <= nT
                for j = 1:blockSize
                    % This assumes a specific mapping, preserving it exactly.
                    if (i * blockSize) + j <= totalSize 
                        FCOPT_TM(rowStart+j-1, (i * blockSize) + j) = WPP_transition_matrix(i, i+1);
                    end
                end
            end
            
            % Assign the left side diagonal values if not in the first block
            if i > 1
                for j = 1:blockSize
                    % This assumes a specific mapping, preserving it exactly.
                    if ((i-2) * blockSize) + j > 0 && ((i-2) * blockSize) + j <= totalSize
                        FCOPT_TM(rowStart+j-1, ((i-2) * blockSize) + j) = WPP_transition_matrix(i, i-1);
                    end
                end
            end
        end
        
        % Find combined probabilities (convolution)
        % Output table initialization
        P_combined = [];
        prob_combined = [];
        
        % Calculate combined states
        for i = 1:size(FCOPT, 1) % FCOPT(i,1) is outage capacity from turbine failures (0 to max_capacity)
            for j = 1:size(WT_States, 1) % WT_States(j,1) is power output from wind (in MW)
                P_new = WT_States(j, 1) * FCOPT(i, 1) / wT; 
                prob_new = FCOPT(i, 2) * WT_States(j, 2); % Combined probability
                P_combined = [P_combined; P_new];
                prob_combined = [prob_combined; prob_new];
            end
        end
        
        % Initialize vectors for storing the sums
        n = size(FCOPT_TM, 1); % Number of rows (and columns) in the matrix
        lambda_plus_combined = zeros(n, 1);
        lambda_minus_combined = zeros(n, 1);
        
        % Loop through each row to compute sums to the right and left of the diagonal
        for i = 1:n
            for j = 1:n
                % P_combined represents a form of "combined outage equivalent power".
                if P_combined(i)<P_combined(j)
                    lambda_plus_combined(i) = lambda_plus_combined(i) + FCOPT_TM(i,j);  % Sum values to the right of the diagonal
                elseif P_combined(i)>P_combined(j)
                    lambda_minus_combined(i) = lambda_minus_combined(i) + FCOPT_TM(i,j);  % Sum values to the left of the diagonal
                end
            end
        end
        
        % Create final combined table
        FCOPT = [P_combined, prob_combined, lambda_plus_combined, lambda_minus_combined];
        
%% Rounding and Final Output Table Generation
% This concluding section takes the detailed combined probability table and rounds it to a user-specified number of states. It then recalculates the probabilities, transition rates, and frequencies for these rounded states, producing the final `rounded_FCOPT` (Farm Capacity Output Probability Table) and `FCOPT_TM_rounded` (Transition Matrix of the Rounded Table).

        % Define the bins for output capacities
        if strcmp(num_rounding_states, 'exact')
            unique_values = unique(P_combined);
            num_rounding_states=length(unique_values); % define exact rounding states
            bins = unique_values'; % Ensure bins is a row vector for linspace replacement
        else
            % Max possible capacity for the bins is total installed capacity (nT*wT)
            bins = linspace(0, 1, num_rounding_states) * wT * nT;
        end
        
        % Initialize the rounded table
        rounded_P = bins';
        rounded_prob = zeros(length(bins), 1);
        rounded_freq_plus = zeros(length(bins), 1);
        rounded_freq_minus = zeros(length(bins), 1);
        indexes_accumulated = cell(length(bins), 1); % Cell array to store indexes for each bin
        
        % Assign probabilities and rates to the closest bins
        for i = 1:size(FCOPT, 1)
            P_current = FCOPT(i, 1);
            
            % Handle the maximum capacity as an exception
            if P_current == wT * nT && ~strcmp(num_rounding_states, 'exact')
                bin_idx = length(bins); % Always add to the last bin for exact matches
            else
                % Find the closest bin, but ensure no rounding upwards to the max capacity
                differences = abs(rounded_P - P_current);
                % if ~strcmp(num_rounding_states, 'exact') % Only apply this exception if not 'exact'
                %     differences(end) = inf; % Set the last bin's difference to infinity unless exact match
                % end
                [~, bin_idx] = min(differences);
            end
            
            % Aggregate probabilities and transition rates to the chosen bin
            rounded_prob(bin_idx) = rounded_prob(bin_idx) + FCOPT(i, 2);
            rounded_freq_plus(bin_idx) = rounded_freq_plus(bin_idx) + FCOPT(i, 2)*FCOPT(i, 3);
            rounded_freq_minus(bin_idx) = rounded_freq_minus(bin_idx) + FCOPT(i, 2)*FCOPT(i, 4);
            indexes_accumulated{bin_idx} = [indexes_accumulated{bin_idx}, i];
        end
        
        % Correction for transition rates insided rounded states
        products_up = cell(length(bins), 1); 
        products_down = cell(length(bins), 1); 
        sum_products_up = zeros(length(bins), 1); 
        sum_products_down = zeros(length(bins), 1); 

        for i = 1:length(bins)
            indexes = indexes_accumulated{i}; % indexes in this bin
            products_up{i} = zeros(length(indexes), length(indexes)); % Initialize matrix of upward products for this bin
            products_down{i} = zeros(length(indexes), length(indexes)); % Initialize matrix of downward products for this bin
            
            for j = 1:length(indexes)
                idx_j = indexes(j);
                for k = 1:length(indexes)
                    idx_k = indexes(k);
                    if idx_j ~= idx_k % Ensure we don't include self-transitions
                        if FCOPT(idx_j, 1) < FCOPT(idx_k, 1)
                            products_up{i}(j, k) = FCOPT(idx_j, 2) * FCOPT_TM(idx_j, idx_k);
                        elseif FCOPT(idx_j, 1) > FCOPT(idx_k, 1)
                            products_down{i}(j, k) = FCOPT(idx_j, 2) * FCOPT_TM(idx_j, idx_k);
                        end
                    end
                end
            end
            
            % Sum all products in the matrix for this bin
            sum_products_up(i) = sum(sum(products_up{i}));
            sum_products_down(i) = sum(sum(products_down{i}));
        end
        
        % Display the products for each bin
        for i = 1:length(bins)
            rounded_freq_plus(i) = rounded_freq_plus(i) - sum_products_up(i);
            rounded_freq_minus(i) = rounded_freq_minus(i) - sum_products_down(i) ;
        end
        
        % Create final rounded table
        % Handling division by zero for rates where rounded_prob is 0.
        temp_rounded_lambda_plus = rounded_freq_plus./rounded_prob;
        temp_rounded_lambda_plus(isnan(temp_rounded_lambda_plus) | isinf(temp_rounded_lambda_plus)) = 0;

        temp_rounded_lambda_minus = rounded_freq_minus./rounded_prob;
        temp_rounded_lambda_minus(isnan(temp_rounded_lambda_minus) | isinf(temp_rounded_lambda_minus)) = 0;

        % Construct rounded_FCOPT exactly as in your script, including flipud
        rounded_FCOPT = [round(rounded_P), rounded_prob, temp_rounded_lambda_plus, temp_rounded_lambda_minus, ((temp_rounded_lambda_plus)+(temp_rounded_lambda_minus)).*rounded_prob];
        rounded_FCOPT=flipud(rounded_FCOPT);
        
        % Create the transition matrix 'FCOPT_TM_rounded'
        freq_rounded = zeros(length(bins), length(bins));
        FCOPT_TM_rounded = zeros(length(bins), length(bins));
        
        % Calculate the transition rates for the rounded states
        for i = 1:length(bins)
            indices_i = indexes_accumulated{i}; % indices for state i
            for j = 1:length(bins)
                if i ~= j
                    indices_j = indexes_accumulated{j}; % indices for state j
                    % Sum transitions from all states in indices_i to all states in indices_j
                    for idx_i = indices_i
                        for idx_j = indices_j
                            freq_rounded(i, j) = freq_rounded(i, j) + FCOPT_TM(idx_i, idx_j)*FCOPT(idx_i, 2); % FCOPT here refers to the intermediate FCOPT before rounding
                        end
                    end
                end
            end
        end
        
        for i = 1:length(bins)
            % Ensure rounded_prob(i) is not zero before division
            if rounded_prob(i) > 0
                for j = 1:length(bins)
                    if i ~= j
                         FCOPT_TM_rounded(i, j) = FCOPT_TM_rounded(i, j) + freq_rounded(i, j)/rounded_prob(i);
                    end
                end
            end
        end
        
        % Mirror the matrix horizontally
        FCOPT_TM_rounded = flip(FCOPT_TM_rounded, 2);
        % Mirror the matrix both horizontally and vertically
        FCOPT_TM_rounded = flip(FCOPT_TM_rounded, 1);

end