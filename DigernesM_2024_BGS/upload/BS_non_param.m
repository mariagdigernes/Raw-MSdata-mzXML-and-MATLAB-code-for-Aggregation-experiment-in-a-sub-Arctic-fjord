%   SUMMARY       : BS_NON_PARAM is a function to perform non-parametric 
%                   tests on two independent samples (Mann-Whitney U test
%                   or Fligner-Policello test).
%                   MATLAB function RANKSUM is used for Mann-Whitney U test
%                   (https://www.mathworks.com/help/stats/ranksum.html).
%                   Fligner-Policello test is a test of two combined random 
%                   variables with continuous cumulative distribution [1]. 
%                   It is a robust rank-order test of treatment effects 
%                   for populations which assumes neither normality nor 
%                   homoscedasticity.
%                   Additonally, to reduce the false-positive rate, 
%                   this function adopts a boot-strapped Monte-Carlo 
%                   method to dynamically approximate the critical 
%                   value of the test statistic when sample sizes < 100
%                   and the test statistic may not follow normal distribution[2].
%	               
%   LIMITATIONS   : 1. This function is only valid for continuous data.
%
%   INPUT         : x,y   = 1D vectors of independent samples
%                   alpha = alpha value
%                   tail  = 'both' / 'right' / 'left'
%                   type  = 1 -> Wilcoxon's ranksum test, 2 -> Fligner-Policello test
%                   p-type = 'approx' or 'exact'
%                 
%   OUTPUT        : p     = p-value of appropriate test of significant
%                           difference (2-tailed)
% 
%   NOTE          : 1. This function requires Statistics and Machine Learning Toolbox
%                   2. Written and tested in MATLAB R2020a
% 
%   EXAMPLE       : (Call from MATLAB terminal)
%	            >> x = normrnd(5,0.5,1000,1);
%                   >> y = normrnd(10,0.5,1000,1);
%	            >> p = RROtest(x,y,0.05);
%                   >> p
%                   ans =
%                   0.0000
%
%   REFERENCES    : 1. Fligner, M.A. and Pollicello, G.E. III. (1981). 
%                      Robust Rank Procedures for the Behrens-Fisher Problem. 
%                      Journal of the American Statistical Association. 
%                      76, 162, 168.
%                  2. Boos, D.D. and Brownie, C. (1988). 
%                     Bootstrap p-Values for Tests of Nonparametric Hypotheses.
%                     Institute of Statistics Mimeo Series No. 1919, 
%                     North Carolina State University.

function p = BS_non_param(x, y, alpha, tail, type, ptype)

    % seed the random number generator after storing starting state
    s = rng;
    rng default;
    
    % PRELIMINARIES   
    % check for proper inputs
    if ~isvector(x) || ~isvector(y) 
        error('each sample (x and y) must be 1 dimensional');
    end

    if ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
        error('alpha value must be a scaler in range [0 1]');
    end
    
    if ~(type == 1 || type == 2)
        error('type should be 1 or 2')
    end
    
    % Check value of 'tail'
    tail = internal.stats.getParamVal(tail, {'both'  'right'  'left'}, '''tail''');
 
    % Check value of 'p-type'
    ptype = internal.stats.getParamVal(ptype, {'approx','exact'}, '''ptype''');
    
    % Remove missing data
    x = x(~isnan(x));
    y = y(~isnan(y));
    if isempty(x) || isempty(y)
        error('Not Enough Data');
    end    
        
    % PART1 (Obtain Mann-Whitney test statistic for data)
    nx = length(x);
    ny = length(y);
    x = x(:);                                                               % ensure columns
    y = y(:);
    if type == 1
        [~,~,stats] = ranksum(x,y);
        if strcmp(ptype,'approx')
            U = stats.zval;
        else
            U = (stats.ranksum - (ny * (ny + 1) / 2)) / (nx * ny);
        end
    else
        U = calc_statisticFP(x,y);
    end
    
    % PART2 (Obtain test statistics for bootstrap data from best 
    %        distribution fits to sample data if sample sizes of 
    %        x and y are <= 100. Otherwise assume standard normal 
    %        distribution of the test statistic)
    if strcmp(ptype,'approx')
        Uvals = random('normal',0,1,1,999);
    else
        all_delta = y - x';
        delta = median(all_delta(:));                                       % Hodges-Lehmann shift operator
        distX = best_distribution(x, alpha);
        distY = best_distribution(y - delta, alpha);

        Uvals = zeros(1,999);
        for i = 1:length(Uvals)
            if length(distX{2}) == 1
                cx = random(distX{1}, distX{2},nx, 1);
            elseif length(distX{2}) == 2
                cx = random(distX{1}, distX{2}(1), distX{2}(2), nx, 1);
            else
                cx = random(distX{1}, distX{2}(1), distX{2}(2), distX{2}(3), nx, 1);
            end

            if length(distY{2}) == 1
                cy = random(distY{1}, distY{2},ny, 1);
            elseif length(distY{2}) == 2
                cy = random(distY{1}, distY{2}(1), distY{2}(2), ny, 1);
            else
                cy = random(distY{1}, distY{2}(1), distY{2}(2), distY{2}(3), ny, 1);
            end       

            % Calculate test statistic
            if type == 1
                [~,~,stats] = ranksum(cx,cy);
                Uvals(i) = (stats.ranksum - (ny * (ny + 1) / 2)) / (nx * ny);
            else
                Uvals(i) = calc_statisticFP(cx,cy);
            end
        end
    end
    
    % Center the test statistic values on the mean
    absdiff = abs(mean(Uvals) - U);
    pleft = (sum(Uvals < (mean(Uvals) - absdiff))) / 1e3;
    pright = (sum(Uvals > (mean(Uvals) + absdiff))) / 1e3;
    if strcmp(tail,'left')
        p = pleft;
    elseif strcmp(tail,'right')
        p = pright;
    else
        p = pleft + pright;
    end
    
    % Restore random number generator initial state
    rng(s);
end

function U = calc_statisticFP(x,y)
    % Compute the placements of each sample element
    U_YXi = zeros(size(x));
    U_XYi = zeros(size(y));   
    
    for i = 1:length(U_YXi)
        U_YXi(i) = sum((x(i) > y));
        U_YXi(i) = U_YXi(i) + 1.5 * sum((x(i) == y));
    end
    
    for i = 1:length(U_XYi)
        U_XYi(i) = sum((y(i) > x));
        U_XYi(i) = U_XYi(i) + 1.5 * sum((y(i) == x));
    end 
    
    %Compute mean placements
    U_YX = mean(U_YXi);
    U_XY = mean(U_XYi);
    
    %Compute a measure of varaince
    Vx = sum((U_YXi - U_YX) .^ 2);
    Vy = sum((U_XYi - U_XY) .^ 2);
    
    U = (length(x) * U_YX - length(y) * U_XY) /...
        (2 * sqrt(Vx + Vy + U_XY * U_YX));
end
