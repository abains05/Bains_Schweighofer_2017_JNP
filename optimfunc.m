function err = optimfunc(fitvar,m_prior_init,sd_prior_init,targets,initang_data)
%function to be minimized when fitting likelihood variance 'fitvar(1)' and prior 
%learning rate 'fitvar(2)' in adaptive Bayesian model from Verstynen and 
%Sabes 2011

%run simulation to generate predicted MAP reach directions from model
theta_map = simulateAdaptiveBayesExp(targets,m_prior_init,sd_prior_init,fitvar(2),fitvar(1),0);

%compare predicted directions to actual data to find error
err = 0;
for i = 1:length(initang_data)
    err = err + circDif(theta_map(i),initang_data(i)) .^2;
end