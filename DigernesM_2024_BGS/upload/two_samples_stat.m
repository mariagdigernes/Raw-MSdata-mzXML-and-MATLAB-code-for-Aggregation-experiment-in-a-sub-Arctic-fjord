%   SUMMARY      : TWO_SAMPLES_STAT is a wrapper to perform unpaired 
%                  two-sample statistical comparison:
%                   1. First tests for parametricity is applied
%                         a. Bartlett's test for homoscedasticity
%                         b. Anderson-Darling test for normality
%                   2. Four cases arise from above:
%                         CASE 1: a = YES, b = YES (Unpaired t-test)
%                         CASE 2: a = NO,  b = YES (Unpaired t-test using
%                                                   Satterthwaite’s approximation)
%                         CASE 3: a = YES, b = NO  (Bootstrapped Mann-Whitney
%                                                   U test)
%                         CASE 4: a = NO,  b = NO  (Bootstrapped
%                                                   Fligner-Pollicello robust
%                                                   rank order test)
%
%   LIMITATIONS  : 1. This function is only valid for continuous data.
%
%   INPUT        : x,y   = 1D vectors of independent samples
%                  alpha = alpha value
%                  tail  = 'both' / 'right' / 'left'
%                 
%   OUTPUT       : p     = p-value of appropriate test of significant
%                          difference (2-tailed)
%
%   NOTE         : 1. The Mann-Whitney U test (which is often used as a
%	              non-parametric alternative to independent t-test) has
%	              a high false positive rate when the two samples have
%	              unequal higher order moments (e.g. unequal variance)
%	              [1]. Fligner and Pollicello's robust rank order test
%	              corrects this short-coming [2].	                  
%                  2. However, when the sample distributions cannot be assumed 
%		      to be normal, they are often fallaciously assumed to be 
%		      at least symmetric which is rarely the case. In fact in 
%		      such cases, the absolute performance of both the above 
%		      non-parametric tests is poor. Additionally, the convergence 
%      		      of the test statistic to standard normal distribution is 
% 		      often slow (as the sample sizes increase) [2,3].
%                  3. Hence, this function adopts a boot-strapped Monte-Carlo 
%                     method to dynamically approximate the critical values 
%                     for the test statistic of the non-parametric tests. 
%                     This technique has the advantage of making no
%                     assumptions about data distributions and therefore
%                     significantly reduce the bias in these tests [3].
%	           4. This function requires Statistics and Machine Learning Toolbox
%                  5. Written and tested in MATLAB R2020a
% 
%	EXAMPLE : (Call from MATLAB terminal)
%	          >> x = normrnd(5,0.5,1000,1);
%                 >> y = normrnd(10,0.5,1000,1);
%                 >> p = two_samples_stat(x,y,0.05);
%                 >> p
%                 ans =
%                 0.0000
%
%  REFERENCES    : 1. Fligner, M.A. and Pollicello, G.E. III. (1981). 
%                     Robust Rank Procedures for the Behrens-Fisher Problem.
%                     Journal of the American Statistical Association. 
%                     76, 162, 68.
%                  2. Feltovich, N. 
%                     Nonparametric Tests of Differences in Medians: 
%                     Comparison of the Wilcoxon-Mann-Whitney and Robust Rank-Order Tests. 
%                     Experimental Economics 6, 273 297 (2003). 
%                     https://doi.org/10.1023/A:1026273319211
%                  3. Boos, D.D. and Brownie, C. (1988). 
%                     Bootstrap p-Values for Tests of Nonparametric Hypotheses. 
%                     Institute of Statistics Mimeo Series No. 1919, 
%                     North Carolina State University.
                       
function p = two_samples_stat(x,y,alpha,tail)
   
    % Bartlett's test for homoscedasticity
    hv = vartest(x,var(y));
    if hv == 1
        flag1 = 1;
    else
        flag1 = 0;
    end
    
    % turn off min p-value warning for AD-test
    warning('off','all');
    
    % Anderson-Darling test for normality
    if length(x) > 4
        h1 = adtest(x);
    else
        h1 = 1;
    end
    
    if length(y) > 4
        h2 = adtest(y);
    else
        h2 = 1;
    end
    
    if h1 == 0 && h2 == 0
        flag2 = 0;
    else
        flag2 = 1;
    end
    warning('on','all');
    
    % Apply appropriate statistical test based on the above verifications
    if flag1 == 0 && flag2 == 0
        [~,p] = ttest2(x,y,'alpha',alpha,'tail',tail);
    end
    
    if flag1 == 1 && flag2 == 0
        [~,p] = ttest2(x,y,'Vartype','unequal','alpha',alpha,'tail',tail);
    end
    
    if flag1 == 0 && flag2 == 1
        p = BS_non_param(x,y,alpha,tail,1);
    end
    
    if flag1 == 1 && flag2 == 1
        p = BS_non_param(x,y,alpha,tail,2);
    end
end
