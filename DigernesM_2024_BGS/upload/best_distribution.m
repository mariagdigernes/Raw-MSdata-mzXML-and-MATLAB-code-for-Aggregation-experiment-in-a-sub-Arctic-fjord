%   SUMMARY      : FIND_BEST_DIST finds best-fitting distribution using MATLAB's function 
%	           MLE ordered by log liklihood. 24 continuous distributions are tested.
%                  Use matlab function RANDOM to sample from the best
%                  obtained distribution for Monte Carlo simulations etc.
%   
%   LIMITATIONS  : 1. This function is only valid for continuous data.
%	
%   INPUT        : data = 1D vector of sample data, alpha = alpha value for mle
%                 
%   OUTPUT       : cell array {'best fit name', [parameter1 , parameter2, parameter3]}
%                  Note that depending on the type of the distribution, the
%                  number of parameters can be 1, 2 or 3
%
%   NOTE         : 1. If X contains negative values, only the Normal distribution is fitted.
%	              Also, if X contains values > 1 the Beta distribution is not fitted. If X
%	              contains 0 some distributions are not fitted.
%	           2. This function requires Statistics and Machine Learning Toolbox
%                  3. Adapted with modifications from Francisco de Castro (2020). fitmethis 
%                     (https://www.mathworks.com/matlabcentral/fileexchange/40167-fitmethis), 
%                     MATLAB Central File Exchange. Retrieved July 18, 2020.
%                  4. Written and tested in MATLAB R2020a
% 
%   EXAMPLE      : (Call from MATLAB terminal)
%	           >> x = normrnd(5,0.5,1000,1);
%	           >> best_dist = find_best_dist(x);
%                  >> best_dist{1}
%                  ans =
%                  'Normal'
%                  >> best_dist{2}
%                  ans = 
%                  4.99 0.494

function best_dist = best_distribution(data, alpha)
    warning off
    % storage variables
    best_dist = {'',[]};
    
    % check for proper inputs
    if ~isvector(data)
        error('data must be 1 dimensional');
    end
    
    if ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
        error('alpha value must be a scaler in range [0 1]');
    end
    
    % Distributions
    dist_list = {'Beta'; 'BirnbaumSaunders'; 'Burr'; 'Exponential'; ...
                 'Extreme Value'; 'Gamma'; 'Generalized Extreme Value';...
                 'Generalized Pareto'; 'Geometric'; 'HalfNormal'; 
                 'InverseGaussian'; 'Loglogistic'; 'Lognormal';'Logistic';...
                 'Uniform'; 'Rayleigh'; 'beta'; 'Nakagami'; 'Normal';...
                 'Rician'; 'Stable'; 'tLocationScale'; 'unif'; 'Weibull'};
    must_be_pos = [2 8 11 12 13 18 20 24];
         
    % Fit only normal distributions for data with negative values
    if min(data) < 0
        phat = mle(data,'distribution','normal','alpha',alpha);
        best_dist(1) = {'normal'};
        best_dist(2) = {phat};
    else      
         LL = -Inf;                                                         %Initialize log liklihood (LL)
         best_dist_idx = -999;                                              %Intialize best_dist index to dummy value
         
        % Fit each distribution in dist_list
        for j= 1:length(dist_list)
            % Check: if values > 1 for beta
            if strcmp('beta',dist_list{j}) && max(data) > 1

            % Check: if values > 0 for some distr
            elseif ismember(j,must_be_pos) && min(data) == 0

            % All other cases
            else
                try
                    phat = mle(data,'distribution',dist_list{j},'alpha',alpha);
                    if length(phat) == 1
                        cpdf = pdf(dist_list{j}, data, phat);
                    elseif length(phat) == 2
                        cpdf = pdf(dist_list{j}, data, phat(1), phat(2));
                    else
                        cpdf = pdf(dist_list{j}, data, phat(1), phat(2), phat(3));
                    end
                    cLL =  sum(log(cpdf(cpdf > 0 & ~isinf(cpdf))));
                    if cLL > LL
                        best_dist_idx = j;
                        LL = cLL;
                        phat_best = phat;
                    end
                catch
                    %This distribution cannot be fitted to the data
                end
            end
        end

        best_dist(1) = dist_list(best_dist_idx);
        best_dist(2) = {phat_best};
    end
end

