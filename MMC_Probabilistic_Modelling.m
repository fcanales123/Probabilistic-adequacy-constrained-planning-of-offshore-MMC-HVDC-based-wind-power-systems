function [n_optimum, T_optimum, min_UAV_value, lambda_arm_optimum, A_mmc_optimum,A_mmc_values] = MMC_Probabilistic_Modelling( ... 
    k_max, V_DC_link, V_IGBT, SF, ...
    pi_v, lambda_igbt, lambda_cap, lambda_diode, ...
    lambda_con, MTTR_con, UAV_0, redundancy, MTTR_arm, MA_Cost, A_mmc_min, ...
    n_max_multiplier, T_values_array)
%MMC_Probabilistic_Modelling Calculates the optimum maintenance interval and
%redundant submodule count for an MMC-based wind power plant.
%
%   [n_optimum, T_optimum, min_UAV_value] = MMC_Probabilistic_Modelling( ...
%   WPP_CAP, k_max, V_DC_link, V_IGBT, SF, ...
%   pi_v, lambda_igbt, lambda_cap, lambda_diode, ...
%   lambda_con, MTTR_con, UAV_0, redundancy, MTTR_arm, MA_Cost, A_mmc_min, ...
%   n_max_multiplier, T_values_array)
%
%   This function determines the optimal number of redundant submodules (n)
%   and the optimal maintenance interval (T) to minimize the Unit Annualized
%   Value (UAV) while ensuring a minimum overall MMC availability.
%
%   Inputs:
%   - k_max (double): Maximum number of submodules that can fail for arm to be operational.
%   - V_DC_link (double): DC link voltage of the MMC.
%   - V_IGBT (double): Voltage rating of a single IGBT.
%   - SF (double): Safety Factor for submodule voltage.
%   - pi_v (double): Operational probability of submodule in a specified state.
%   - lambda_igbt (double): Failure rate of an IGBT.
%   - lambda_cap (double): Failure rate of a capacitor.
%   - lambda_diode (double): Failure rate of a diode.
%   - lambda_con (double): Failure rate at the converter level.
%   - MTTR_con (double): Mean Time To Repair for converter level failures.
%   - UAV_0 (double): Annualized cost of a single submodule (SM_Cost).
%   - redundancy (char array): Type of redundancy ('SBY' for Standby, 'KON' for k-out-of-n).
%   - MTTR_arm (double): Mean Time To Repair for an arm failure.
%   - MA_Cost (double): Cost of a single maintenance action.
%   - A_mmc_min (double): Minimum acceptable availability for the overall MMC.
%   - n_max_multiplier (double, optional): Multiplier for calculating maximum redundant submodules (default: 0.1).
%   - T_values_array (vector, optional): Array of maintenance intervals to evaluate in days (default: [0.5, ..., 10]).
%
%   Outputs:
%   - n_optimum (double): Optimal number of redundant submodules.
%   - T_optimum (double): Optimal maintenance interval in days.
%   - min_UAV_value (double): Minimum Unit Annualized Value achieved.

    % Set default values for optional inputs if not provided
    if nargin < 17
        n_max_multiplier = 0.1;
    end
    if nargin < 18
        T_values_array = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10];
    end
    
    %% Reliability Function Definitions
    % These nested functions define the reliability and expected uptime
    % for different redundancy schemes of the MMC arms.
    
    % Reliability function for k-out-of-n redundancy
    % R_X_k_over_n(T, l, lambda_SM, n) calculates the probability that at least
    % 'l' out of 'n+l' submodules are operational over a period 'T'.
    R_X_k_over_n = @(T, l, lambda_SM, n) sum(arrayfun(@(i) nchoosek(n+l, i) * exp(-i * lambda_SM * T) * (1 - exp(-lambda_SM * T))^(n+l-i), l:(n+l)));
    
    % Integral of R_X_k_over_n(T) from 0 to T
    % This integral is used in the calculation of Mean Time To Failure (MTTF).
    integral_R_X_k_over_n = @(T, l, lambda_SM, n) sum(arrayfun(@(i) ...
        (1 / (i * lambda_SM)) * (1 - sum(arrayfun(@(k_idx) ... % Renamed 'k' to 'k_idx' to avoid conflict with loop variable 'k' in outer scope
        nchoosek(n + l, i + k_idx) * exp(-(i + k_idx) * lambda_SM * T) * ...
        (1 - exp(-lambda_SM * T))^(n + l - (i + k_idx)), 0:(n + l - i)))), l:(n + l)));
    
    % Expected uptime for k-out-of-n system before the first failure requiring repair
    E_X_k_over_n_star = @(T, l, lambda_SM, n) ...
        integral_R_X_k_over_n(T, l, lambda_SM, n) / (1 - R_X_k_over_n(T, l, lambda_SM, n));
    
    % Reliability function for Standby redundancy
    % R_X_standby(T, l, lambda_SM, n) calculates the probability that a standby
    % system with 'l' active and 'n' standby submodules remains operational over 'T'.
    R_X_standby = @(T, l, lambda_SM, n) exp(-l * lambda_SM * T) * sum(arrayfun(@(k_idx) ((l * lambda_SM * T)^k_idx / factorial(k_idx)), 0:n)); % Renamed 'k' to 'k_idx'
    
    % Integral of R_X_standby(T) from 0 to T
    % This integral is used in the calculation of Mean Time To Failure (MTTF).
    integral_R_X_standby = @(T, l, lambda_SM, n) sum(arrayfun(@(k_idx) (1 / (l * lambda_SM)) - ... % Renamed 'k' to 'k_idx'
        sum(arrayfun(@(p) ((l * lambda_SM * T)^p * exp(-l * lambda_SM * T)) / (factorial(p) * l * lambda_SM), 0:k_idx)), 0:n));
    
    % Expected uptime for standby system before the first failure requiring repair
    E_X_standby_star = @(T, l, lambda_SM, n) ...
        integral_R_X_standby(T, l, lambda_SM, n) / (1 - R_X_standby(T, l, lambda_SM, n));
    
    %% MMC Parameter Calculations
    % Calculate critical parameters for the MMC, such as the minimum required
    % submodules and overall submodule failure rates.
    
    % Least necessary number of submodules for an arm to function (based on voltage requirements)
    l = ceil(k_max * V_DC_link / (V_IGBT * SF));
    
    % Total submodule failure rate based on component failure rates
    lambda_SM_values = pi_v * (2 * lambda_igbt + lambda_cap) + lambda_diode;
    
    % Converter-level availability (accounting for balance and control system failures)
    A_con = 1 / (1 + lambda_con * MTTR_con);
    
    % Range of values for 'l' (least necessary submodules) - kept for consistency, but here it's a single value
    l_values = l;
    
    % Maximum number of redundant submodules to consider, expressed as a percentage of 'l'
    n_max = ceil(n_max_multiplier * l_values);
    
    % Cost of a single submodule (Unit Annualized Value of the submodule)
    SM_Cost = UAV_0;
    
    % Values for n (number of redundant submodules) to iterate through
    n_values = linspace(0, n_max, n_max + 1);
    
    % Input T range (Maintenance interval in days)
    T_values = T_values_array;
    
    %% Optimization Loop
    % Iterates through all combinations of 'n' (redundant submodules) and
    % 'T' (maintenance interval) to find the configuration that minimizes
    % the UAV while meeting the minimum availability requirement.
    
    % Preallocate arrays for storing results for all n and T combinations
    A_mmc_values = zeros(length(n_values), length(T_values));
    UAV_values = zeros(length(n_values), length(T_values));
    
    % Variables to store the optimal values
    min_UAV_value = Inf; % Initialize with a very high value
    n_optimum = NaN;     % Initialize as Not-a-Number
    T_optimum = NaN;     % Initialize as Not-a-Number
    A_mmc_optimum = NaN; % Store the availability at the optimum point
    lambda_arm_optimum = NaN; % Store the arm failure rate at the optimum point
    
    % Loop through each value of n (number of redundant submodules)
    for i = 1:length(n_values)
        % Loop through each value of T (maintenance interval)
        for k = 1:length(T_values)
            % Assign current iteration's parameters
            lambda_SM = lambda_SM_values; % Submodule failure rate
            l = l_values;                 % Least necessary submodules
            n = n_values(i);              % Current number of redundant submodules
            T = T_values(k);              % Current maintenance interval
            
            % Calculate Mean Time To Failure (MTTF) for an arm based on redundancy type
            if strcmp(redundancy, 'SBY')
                MTTF_arm = E_X_standby_star(T, l, lambda_SM, n);
            elseif strcmp(redundancy, 'KON')
                MTTF_arm = E_X_k_over_n_star(T, l, lambda_SM, n);
            else
                error('Invalid redundancy type. Please choose ''SBY'' or ''KON''.');
            end
            
            % Calculate reliability parameters for a single MMC arm
            a = 1 / MTTF_arm; % Failure rate of an arm
            b = 1 / MTTR_arm; % Repair rate of an arm
            A_arm = b / (a + b); % Availability of a single arm
            
            % Calculate overall MMC availability (assuming 6 arms in the MMC)
            A_mmc = A_arm^6 * A_con;
            A_mmc_values(i, k) = A_mmc; % Store for analysis if needed
            
            % Calculate Unit Annualized Value for the current configuration
            % This includes the cost of redundant submodules and maintenance costs.
            UAV_values(i, k) = n * SM_Cost + MA_Cost / T;
           
            % Check if the current configuration meets the minimum required
            % MMC availability AND has a lower UAV than the current minimum.
            if A_mmc >= A_mmc_min && UAV_values(i, k) < min_UAV_value
                min_UAV_value = UAV_values(i, k); % Update minimum UAV
                n_optimum = n;                    % Update optimal n
                T_optimum = T;                    % Update optimal T
                A_mmc_optimum = A_mmc;            % Store optimal A_mmc (for internal debug/info)
                lambda_arm_optimum = a;           % Store optimal arm failure rate (for internal debug/info)
            end
        end
    end

    % Display warnings if no optimal solution was found
    if isnan(n_optimum)
        warning('MMC_Probabilistic_Modelling:NoOptimalSolution', ...
                'No (n, T) combination found that satisfies the minimum availability (A_mmc_min) and minimizes UAV.');
        fprintf('Consider adjusting input parameters, especially A_mmc_min, or extending the search ranges for n and T.\n');
    end

end