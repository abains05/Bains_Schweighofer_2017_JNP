function [sd_likelihood,B,fval,exitflag] = fitBayesAdaptive(m_prior_init,sd_prior_init,targets,initang_data,algo)
%uses initial reach angle data (specified in initang_data), as well as the 
%initial values of the prior mean 'm_prior_init' and prior SD 'sd_prior_init', and 
%actual targets to fit the likelihood SD and learning rate in the adaptive 
%Bayesian model from Vertsynen and Sabes 2011. 'algo' determines which
%optimization algorithm to use in fmincon.

%initialize params to be fit
B_init = 0.8;
sd_likelihood_init = 12;

%find fit using fmincon
f = @(fitvar)optimfunc(fitvar,m_prior_init,sd_prior_init,targets,initang_data);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-6,'Algorithm',algo);
[fitvar,fval,exitflag] = fmincon(f,[sd_likelihood_init; B_init],[],[],[],[],[0; 0],[inf; 1],[],options);
sd_likelihood = fitvar(1);
B = fitvar(2);
