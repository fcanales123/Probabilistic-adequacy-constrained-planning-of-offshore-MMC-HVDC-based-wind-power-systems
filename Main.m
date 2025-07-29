clear;
tic;             % Start the timer
set(findall(gcf,'type','axes'),'fontsize',12,'FontName', 'Times New Roman')
set(findall(gcf,'type','text'),'fontSize',12,'FontName', 'Times New Roman')

turbine='IEA'; % Turbine model used 'IEA' or 'DTU'
redundancy='SBY'; % Redundancy strategy used 'SBY'-stand-by or 'KON' k-out-of-n
load_model= 'comulative'; %Load model used 'exact' or 'comulative'

% Wind Power Capacities for loop evaluation
if strcmp(turbine,'IEA')
    WPP_CAP_loop = 600:15:2010; %naive=600, Optima_5states=825 , Optimal_exact=960 
elseif strcmp(turbine,'DTU')
    WPP_CAP_loop = 600;
end

% MMC availability constraint for loop evaluation
A_mmc_min_loop = linspace(0.98,0.9991,length(WPP_CAP_loop));  % Minimum required availability of MMC (max 0.99913675), naive=0.995, , Optimal_5states=0.984267021276596, Optimal_exact=0.984267021276596

ConvUSD = 1.3; % Conversion ratio UK pounds to USD
offshore_distance = 100; % Distance offshore in km (a number from 10-100)
num_rounding_states = 5; % Number of states rounded to in the WF model use 'exact' if no rounding required
num_wind_states=5; % Number of states rounded to in the WT (wind speed) model 

% MMC Parameters
V_IGBT=6.5; % KV switch voltage
k_max=1.1; SF=0.6; % k_max voltage ripple factor and SF safety factor

lambda_igbt = 8.76e-4; lambda_cap = 1.752e-3; lambda_diode = 6.123e-3; % SM components failure rates in [occ/year]
pi_v = 1; % failure rate voltage factor only applies to k-out-of-n leave 1 for stand-by

lambda_con= 2.628e-2; MTTR_con=12/365; % Converter level failures in [occ/year]
MTTR_arm = 12 / 365; % Mean time to repair for arm in [occ/year]

UAV_0=2.975; % Uniform anual value per SM with no maintenance in [k$]
MA_Cost=60; % Cost of a maintenance activity in [k$]
n_max_multiplier=0.1; % Maximum percentage redundant SM with respect to required active SMs
T_values_array = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]; % Array of maintenance intervals to evaluate in years. 

% DCTCP Block Parameters
V_DC_link=525; % KV link in DC 
lambda_ACB=0.025; lambda_CT=0.037; lambda_CPD=0.088; lambda_ACF=0.2; lambda_CR=0.116; lambda_DCB=0.121; % (occ/y) 
lambda_DCL=0.071; % (occ/y*100km)

MTTR_ACB=168; MTTR_CT=1580; MTTR_CPD=6.05; MTTR_CR=256; MTTR_ACF=6;  % (h)
MTTR_DCL=1440; % (h/100km)
mu_DCB=121.7; % (occ/y*100km)

cable_length=100; % (km) Total from receiving to sending end
nSEDCTCPs= 2; % integer number for number of sending ends in MTDC
nREDCTCPs= 2; % integer number for number of receiving ends in MTDC

e = 0.5; % Exposure factor for frequency-duration (exact) model
Annual_peak=2850; % Annual peak load for Load duration curve

% Pre-allocate the output arrays before the parfor loop
EWSUP = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));
EENS = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));
ESWE = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));
Cost = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));
Sustainability = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));
Reliability = zeros(length(WPP_CAP_loop), length(A_mmc_min_loop));

parfor x=1:length(WPP_CAP_loop)
    % Initialize WPP_OM_COST here for each parfor iteration
    WPP_OM_COST = 0; % Or a more appropriate initial value
    
    % Create temporary local variables for results within this parfor iteration's slice
    % This is good practice for parfor to make sure results are accumulated correctly
    % before assigning them to the sliced output variables.
    local_EWSUP_row = zeros(1, length(A_mmc_min_loop));
    local_EENS_row = zeros(1, length(A_mmc_min_loop));
    local_ESWE_row = zeros(1, length(A_mmc_min_loop));
    local_Cost_row = zeros(1, length(A_mmc_min_loop));
    local_Sustainability_row = zeros(1, length(A_mmc_min_loop));
    local_Reliability_row = zeros(1, length(A_mmc_min_loop));

    for y=1:length(A_mmc_min_loop)
        % WPP Parameters
        WPP_CAP= WPP_CAP_loop(x); % Total connected wind farm power in MW [600, 800||810, 1000||1005, 1200, 1400||1410, 1600||1605, 1800, 2000||2010] DTU||IEA
        [rounded_FCOPT, FCOPT_TM_rounded, WT_States] = Wind_Farm_Probabilistic_Modelling(turbine, WPP_CAP, offshore_distance, num_wind_states, num_rounding_states);   
        
        % % Display results
        % disp('Rounded FCOPT Combined Capacity Table:');
        % disp('Cap. X (MW)   Prob. p(X)   λ+ (occ/day)   λ- (occ/day)   freq. (occ/day)');
        % disp(rounded_FCOPT);
        % 
        % % Display the rounded transition matrix
        % disp('Rounded Transition Matrix (FCOPT_TM_rounded):');
        % disp(FCOPT_TM_rounded);
        
        % % Evaluating limiting state probabilities FCOPT_TM_rounded for verification
        % Delta_t= 1; % 1 y (doesn't matter, leave it 1)
        % % create STP
        % STP_FCOPT_TM_rounded=FCOPT_TM_rounded*Delta_t;
        % 
        % % Adjust the diagonal entries
        % for i = 1:size(STP_FCOPT_TM_rounded, 1)
        %     STP_FCOPT_TM_rounded(i, i) = 1 - sum(STP_FCOPT_TM_rounded(i, :));
        % end
        % 
        % % Calculate the alpha values
        % [V, D] = eig(STP_FCOPT_TM_rounded.');
        % 
        % % Find the index of the eigenvalue that is 1
        % [~, idx] = min(abs(diag(D) - 1));
        % 
        % % Extract the corresponding eigenvector
        % alpha_FCOPT = V(:, idx);
        % % Normalize the eigenvector so that its elements sum to 1
        % alpha_FCOPT = alpha_FCOPT / sum(alpha_FCOPT);
        % 
        % % Display the solution
        % disp('Solution for alpha_FCOPT:');
        % disp(alpha_FCOPT');
        
        % MMC reliability model
        A_mmc_min = A_mmc_min_loop(y); % Minimum required availability of MMC (min 0.9800 max 0.9978)

        [n_optimum, T_optimum, min_UAV_value, lambda_arm_optimum, A_mmc_optimum] = ...
        MMC_Probabilistic_Modelling(k_max, V_DC_link, V_IGBT, SF, ...
        pi_v, lambda_igbt, lambda_cap, lambda_diode, ...
        lambda_con, MTTR_con, UAV_0, redundancy, MTTR_arm, MA_Cost, A_mmc_min, ...
        n_max_multiplier, T_values_array);

        % Display the optimum n and T values
        % fprintf('Optimum n: %.2f\n', n_optimum);
        % fprintf('Optimum T: %.2f\n', T_optimum);
        % fprintf('Optimum SM UAV: %.2f\n', min_UAV_value);

        MTDC_CAP= WPP_CAP; % Use this to round transmission to a different value than WPP_CAP
        
        [SEDCTCPCOPT, REDCTCPCOPT, DCTCP_TM, DCTCPCOPT] = DC_Station_Probabilistic_Modelling( ...
        MTTR_ACB, MTTR_CT, MTTR_CPD, MTTR_ACF, MTTR_CR, MTTR_DCL, ...
        lambda_ACB, lambda_CT, lambda_CPD, lambda_ACF, lambda_CR, ...
        lambda_arm_optimum, lambda_con, A_mmc_optimum, ...
        lambda_DCB, mu_DCB, lambda_DCL, cable_length, ...
        WPP_CAP, MTDC_CAP, nSEDCTCPs, nREDCTCPs);

        % Display the result
        % disp('Rounded sending DCTCPCOPT Combined Capacity Table:');
        % disp('Cap. X (MW)   Prob. p(X)   λ+ (occ/day)   λ- (occ/day)   freq. (occ/day)');
        % disp(SEDCTCPCOPT);
        % 
        % % Display the result
        % disp('Rounded receiving DCTCPCOPT Combined Capacity Table:');
        % disp('Cap. X (MW)   Prob. p(X)   λ+ (occ/day)   λ- (occ/day)   freq. (occ/day)');
        % disp(REDCTCPCOPT);
        
        % % Evaluating limiting state probabilities of DCTCP_TM for verification 
        % Delta_t= 1; % 1 y (doesn't matter, leave it 1)
        % % create STP
        % STP_DCTCP_TM=DCTCP_TM*Delta_t;
        % 
        % % Adjust the diagonal entries
        % for i = 1:size(STP_DCTCP_TM, 1)
        %     STP_DCTCP_TM(i, i) = 1 - sum(STP_DCTCP_TM(i, :));
        % end
        % 
        % % Calculate the alpha values
        % [V, D] = eig(STP_DCTCP_TM.');
        % 
        % % Find the index of the eigenvalue that is 1
        % [~, idx] = min(abs(diag(D) - 1));
        % 
        % % Extract the corresponding eigenvector
        % alpha_DCTCP = V(:, idx);
        % 
        % % Normalize the eigenvector so that its elements sum to 1
        % alpha_DCTCP = alpha_DCTCP / sum(alpha_DCTCP);
        % 
        % % Display the solution
        % disp('Solution for alpha_DCTCP:');
        % disp(alpha_DCTCP');
   
        [MTDCCOPT_merged, MTDC_TM_merged] = MTDC_Probabilistic_Modelling(nSEDCTCPs, nREDCTCPs, DCTCPCOPT, DCTCP_TM, MTDC_CAP);

        % % Display the results
        % disp('Merged Capacity States Table (MTDCCOPT):');
        % disp('Cap. X (MW)   Prob. p(X)   λ+   λ-   freq.');
        % disp(MTDCCOPT_merged);

        % % Evaluating limiting state probabilities of MTDC_TM_merged for verification  
        % 
        % Delta_t= 1; % 1 y (doesn't matter, leave it 1)
        % 
        % % create STP
        % STP_MTDC_TM_merged=MTDC_TM_merged*Delta_t;
        % 
        % % Adjust the diagonal entries
        % for i = 1:size(STP_MTDC_TM_merged, 1)
        %     STP_MTDC_TM_merged(i, i) = 1 - sum(STP_MTDC_TM_merged(i, :));
        % end
        % 
        % % Calculate the alpha values
        % [V, D] = eig(STP_MTDC_TM_merged.');
        % 
        % % Find the index of the eigenvalue that is 1
        % [~, idx] = min(abs(diag(D) - 1));
        % 
        % % Extract the corresponding eigenvector
        % alpha_MTDC_TM_merged = V(:, idx);
        % 
        % % Normalize the eigenvector so that its elements sum to 1
        % alpha_MTDC_TM_merged = alpha_MTDC_TM_merged / sum(alpha_MTDC_TM_merged);
        % 
        % % Display the solution
        % disp('Solution for alpha_MTDC_TM_merged:');
        % disp(alpha_MTDC_TM_merged');
        % 
        [CONVOLVED_SCOPT_merged, CONVOLVED_SCOPT_TM_merged] = Convolution_WPP_MTDC(FCOPT_TM_rounded, MTDC_TM_merged, rounded_FCOPT, MTDCCOPT_merged);
        
        % % Evaluating limiting state probabilities of CONVOLVED_SCOPT_TM_merged for verification  
        % 
        % Delta_t= 1; % 1 y (doesn't matter, leave it 1)
        % 
        % % create STP
        % STP_CONVOLVED_SCOPT_TM_merged=CONVOLVED_SCOPT_TM_merged*Delta_t;
        % 
        % % Adjust the diagonal entries
        % for i = 1:size(STP_CONVOLVED_SCOPT_TM_merged, 1)
        %     STP_CONVOLVED_SCOPT_TM_merged(i, i) = 1 - sum(STP_CONVOLVED_SCOPT_TM_merged(i, :));
        % end
        % 
        % % Calculate the alpha values
        % [V, D] = eig(STP_CONVOLVED_SCOPT_TM_merged.');
        % 
        % % Find the index of the eigenvalue that is 1
        % [~, idx] = min(abs(diag(D) - 1));
        % 
        % % Extract the corresponding eigenvector
        % alpha_CONVOLVED_SCOPT_TM_merged = V(:, idx);
        % 
        % % Normalize the eigenvector so that its elements sum to 1
        % alpha_CONVOLVED_SCOPT_TM_merged = alpha_CONVOLVED_SCOPT_TM_merged / sum(alpha_CONVOLVED_SCOPT_TM_merged);
        % 
        % % Display the results
        % disp('Merged  WPP + MTDC System (CONVOLVED_SCOPT_merged) Capacity Table:');
        % disp('Cap. X (MW)   Prob. p(X)   λ+   λ-   freq.');
        % disp(CONVOLVED_SCOPT_merged);
        % 
        % % Display the solution
        % disp('Solution for alpha_CONVOLVED_SCOPT_TM_merged:');
        % disp(alpha_CONVOLVED_SCOPT_TM_merged');

        units = [
                struct('Capacity', 12, 'States', [12, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.008163265, 0], 'RepairRate', [0, 0.4]),
                struct('Capacity', 12, 'States', [12, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.008163265, 0], 'RepairRate', [0, 0.4]),
                struct('Capacity', 12, 'States', [12, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.008163265, 0], 'RepairRate', [0, 0.4]),
                struct('Capacity', 12, 'States', [12, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.008163265, 0], 'RepairRate', [0, 0.4]),
                struct('Capacity', 12, 'States', [12, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.008163265, 0], 'RepairRate', [0, 0.4]),
                struct('Capacity', 20, 'States', [20, 0], 'Probabilities', [0.9, 0.1], 'FailureRate', [0.053333333, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 20, 'States', [20, 0], 'Probabilities', [0.9, 0.1], 'FailureRate', [0.053333333, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 20, 'States', [20, 0], 'Probabilities', [0.9, 0.1], 'FailureRate', [0.053333333, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 20, 'States', [20, 0], 'Probabilities', [0.9, 0.1], 'FailureRate', [0.053333333, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 50, 'States', [50, 0], 'Probabilities', [0.99, 0.01], 'FailureRate', [0.012121212, 0], 'RepairRate', [0, 1.2]),
                struct('Capacity', 76, 'States', [76, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.012244898, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 76, 'States', [76, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.012244898, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 76, 'States', [76, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.012244898, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 76, 'States', [76, 0], 'Probabilities', [0.98, 0.02], 'FailureRate', [0.012244898, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 100, 'States', [100, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.02, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 100, 'States', [100, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.02, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 100, 'States', [100, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.02, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 155, 'States', [155, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.025, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 155, 'States', [155, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.025, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 155, 'States', [155, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.025, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 155, 'States', [155, 0], 'Probabilities', [0.96, 0.04], 'FailureRate', [0.025, 0], 'RepairRate', [0, 0.6]),
                struct('Capacity', 197, 'States', [197, 0], 'Probabilities', [0.95, 0.05], 'FailureRate', [0.025263158, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 197, 'States', [197, 0], 'Probabilities', [0.95, 0.05], 'FailureRate', [0.025263158, 0], 'RepairRate', [0, 0.48]),
                struct('Capacity', 197, 'States', [197, 0], 'Probabilities', [0.95, 0.05], 'FailureRate', [0.025263158, 0], 'RepairRate', [0, 0.48]),
                % struct('Capacity', 350, 'States', [350, 0], 'Probabilities', [0.92, 0.08], 'FailureRate', [0.020869565, 0], 'RepairRate', [0, 0.24]), 
                struct('Capacity', 400, 'States', [400, 0], 'Probabilities', [0.88, 0.12], 'FailureRate', [0.021818182, 0], 'RepairRate', [0, 0.16]),
                struct('Capacity', 400, 'States', [400, 0], 'Probabilities', [0.88, 0.12], 'FailureRate', [0.021818182, 0], 'RepairRate', [0, 0.16]),
                struct('Capacity', max(CONVOLVED_SCOPT_merged(:,1)), 'States', CONVOLVED_SCOPT_merged(:,1)', 'Probabilities', CONVOLVED_SCOPT_merged(:,2)', 'FailureRate', CONVOLVED_SCOPT_merged(:,4)', 'RepairRate', CONVOLVED_SCOPT_merged(:,3)')
            ];
        % Decomissioned 350 MW of coal fired generation

        % Assume max_capacity and initialization are done similarly
        %Initialization for array sizes
        max_capacity = sum(arrayfun(@(x) x.Capacity, units)); % Define the maximum capacity to initialize arrays
        com_cap = 0; %Initialize an accumulating capacity variable for each i unit added
        
        % Initialize the COPT array and transition rates
        SCOPT = []; % Final System Outage probability table from 0 out to max_capacity (all units out) and all its possible intermediate combinations
        COPT = zeros(max_capacity + 1, 1);% Recursive outage probability table from 0 out to max_capacity 
        % flp_COPT = zeros(max_capacity + 1, 1); % Flipped recursive outage probability table from 0 out to max_capacity
        lambda_plus = zeros(max_capacity + 1, 1);  % Upward transitions
        lambda_minus = zeros(max_capacity + 1, 1);  % Downward transitions
        
        % Start with the system fully available
        COPT(1) = 1;
        
        % Update loop to include derated states
        for i = 1:length(units)
            unit = units(i);
            % Temporary arrays for probabilities and rates
            temp_COPT = zeros(max_capacity + 1, 1);
            temp_lambda_plus = zeros(max_capacity + 1, 1);
            temp_lambda_minus = zeros(max_capacity + 1, 1);
        
            for state_idx = 1:length(unit.States)
                C_i = unit.States(state_idx);
                p_i = unit.Probabilities(state_idx);
                mu_i = unit.RepairRate(state_idx);
                lambda_i = unit.FailureRate(state_idx);
        
                for j = 0:max_capacity
                    if (j + C_i <= max_capacity) && (COPT(j + 1) > 0)
                        % Update probability p(X)
                        temp_COPT(j + C_i + 1) = temp_COPT(j + C_i + 1) + COPT(j + 1) * p_i;
        
                        % Update λ+(X) and λ-(X) sums, actual division to occur after loop
                        temp_lambda_plus(j + C_i + 1) = temp_lambda_plus(j + C_i + 1) + COPT(j + 1) * p_i * (lambda_plus(j + 1) + mu_i);
                        temp_lambda_minus(j + C_i + 1) = temp_lambda_minus(j + C_i + 1) + COPT(j + 1) * p_i * (lambda_minus(j + 1) + lambda_i);
                    end
                end
            end
        
            % Normalize λ+(X) and λ-(X)
            valid_indices = find(temp_COPT > 0);
            lambda_plus(valid_indices) = temp_lambda_plus(valid_indices) ./ temp_COPT(valid_indices);
            lambda_minus(valid_indices) = temp_lambda_minus(valid_indices) ./ temp_COPT(valid_indices);
        
            % Update COPT
            COPT = temp_COPT;
        end
        
        % Write system COPT only in posible states of outage from max_capacity (all units out) to 0 out and all its possible intermediate combinations
        for X = 0:max_capacity
            if COPT(max_capacity - X + 1) ~= 0
                SCOPT = [SCOPT;max_capacity-X,COPT(max_capacity - X + 1),lambda_plus(max_capacity - X + 1),lambda_minus(max_capacity - X + 1),COPT(max_capacity - X + 1)*(lambda_plus(max_capacity - X + 1)+lambda_minus(max_capacity - X + 1))]; 
            end
        end

        SCOPT_trimmed = SCOPT;
        n = size(SCOPT_trimmed, 1);               % Number of rows
        cumulative = 1;                           % Starting cumulative probability
        cumulative_probs = zeros(n, 1);           % Preallocate result column
        threshold = 1e-8;                         % Stop condition
        
        for i = 1:n
            cumulative_probs(i) = cumulative;     % Store current cumulative value
            decrement = SCOPT_trimmed(i,2);       % Value to subtract
            cumulative = cumulative - decrement;
        
            if cumulative < threshold
                % Stop once below threshold
                cumulative_probs = cumulative_probs(1:i);        % Trim result
                SCOPT_trimmed = SCOPT_trimmed(1:i, :);           % Trim original matrix
                break;
            end
        end

        [LOLP, LOLE, LOEE, EIR, Wind_Eshare, EENS_0, EENS_W] = Reliability_Evaluation(load_model, SCOPT_trimmed, Annual_peak, CONVOLVED_SCOPT_merged, SCOPT);
        
        % Assign to local row variables first
        local_EWSUP_row(y) = EENS_0 + EENS_W;
        local_EENS_row(y) = LOEE;
        local_ESWE_row(y) = sum(rounded_FCOPT(:,1) .* rounded_FCOPT(:,2))*24*365 -  (EENS_0 + EENS_W);
        
        if strcmp(turbine,'IEA')
            if offshore_distance < 50 && offshore_distance >= 10    
                WPP_OM_COST = round((EENS_0 + EENS_W)*ConvUSD*(16.98+0.048*(offshore_distance-10)));  
            elseif offshore_distance >= 50
                WPP_OM_COST = round((EENS_0 + EENS_W)*ConvUSD*(18.9+0.437*(offshore_distance-50)));
            end
        elseif strcmp(turbine,'DTU')
            if offshore_distance < 50 && offshore_distance >= 10    
                WPP_OM_COST = round((EENS_0 + EENS_W)*ConvUSD*(29.53+0.065*(offshore_distance-10)));  
            elseif offshore_distance >= 50
                WPP_OM_COST = round((EENS_0 + EENS_W)*ConvUSD*(32.13+0.552*(offshore_distance-50))); 
            end
        end
        local_Cost_row(y) = (min_UAV_value*n_optimum*6*2*(nREDCTCPs+nSEDCTCPs)*1000+WPP_OM_COST)/(EENS_0 + EENS_W); % 6 arms and 2 MMCs
        local_Sustainability_row(y) = Wind_Eshare; 
        local_Reliability_row(y) = EIR*100; 
    end
    % Assign the accumulated row to the sliced output variables
    EWSUP(x,:) = local_EWSUP_row;
    EENS(x,:) = local_EENS_row;
    ESWE(x,:) = local_ESWE_row;
    Cost(x,:) = local_Cost_row;
    Sustainability(x,:) = local_Sustainability_row;
    Reliability(x,:) = local_Reliability_row;
end

rounded_matrix = round(Reliability, 3);
[rows, cols] = find(rounded_matrix == 99.984); % Original reliability evaluated using comulative method 99.9836866363668 and exact method 99.9891731571509

numMatches = length(rows);
selected_costs = zeros(numMatches, 1);

for i = 1:numMatches
    selected_costs(i) = Cost(rows(i), cols(i));
end

% Step 3: Find the index (among selected) with the minimum cost
[~, minIdx] = min(selected_costs);

% Step 4: Get the row and column in the original matrix
minRow = rows(minIdx);
minCol = cols(minIdx);

% Display the result
fprintf('Minimum cost among selected indices is at (%d, %d)\n', minRow, minCol);
fprintf('Cost [$/MWh] = %.4f\n', Cost(minRow, minCol));
fprintf('Wind Energy Penetration [%%] = %.4f\n', Sustainability(minRow, minCol));
fprintf('WPP capacity [MW] = %.4f\n', WPP_CAP_loop(minRow));
fprintf('MMC Availability Constraint [%%] = %.4f\n', A_mmc_min_loop(minCol)*100);

elapsedTime = toc;  % Stop the timer and get elapsed time
fprintf('Elapsed time: %.4f seconds\n', elapsedTime);

% Create the figure
figure;
surf(Sustainability, Cost, Reliability, 'EdgeColor', 'none');

% Set axis labels
xlabel('Sustainability (Wind Energy Penetration [%])');
ylabel('Cost (Offshore O&M [$/MWh])');
zlabel('Reliability (Energy Index EIR [%])');

% % Add colorbar
% colorbar;
colormap(jet);

xlim([min(Sustainability(:)) max(Sustainability(:))]);
ylim([min(Cost(:)) max(Cost(:))]);
zlim([min(Reliability(:)) max(Reliability(:))]);
