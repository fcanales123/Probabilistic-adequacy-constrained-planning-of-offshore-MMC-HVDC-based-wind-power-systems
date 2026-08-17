function copulaModel = Fit_TwoSite_Wind_Copula(turbine, file1, file2)

    %% Load both wind data files
    S1 = load(file1);
    S2 = load(file2);

    f1 = fieldnames(S1);
    f2 = fieldnames(S2);

    D1 = S1.(f1{1});
    D2 = S2.(f2{1});

    %% Select wind-speed column
    if strcmp(turbine,'IEA')
        windCol = 2;
    elseif strcmp(turbine,'DTU')
        windCol = 1;
    else
        error('Invalid turbine type.');
    end

    %% Ensure paired samples
    N = min(size(D1,1), size(D2,1));

    w1 = D1(1:N,windCol);
    w2 = D2(1:N,windCol);

    valid = isfinite(w1) & isfinite(w2);

    w1 = w1(valid);
    w2 = w2(valid);

    N = length(w1);

    %% Convert observations to pseudo-uniform variables
    U1 = (tiedrank(w1) - 0.5)/N;
    U2 = (tiedrank(w2) - 0.5)/N;

    U = [U1 U2];

    %% Fit Gaussian copula
    Rho = copulafit('Gaussian', U);

    %% Store model
    copulaModel.Rho = Rho;
    copulaModel.w1 = w1;
    copulaModel.w2 = w2;

    fprintf('Gaussian copula correlation = %.4f\n', Rho(1,2));

end