function [theta_map,m_prior,sd_prior] = simulateAdaptiveBayesExp(theta,m_prior_init,sd_prior_init,B,sd_likelihood,sensory_uncertainty_flag)
%Given actual targets 'theta', prior learning rate 'B', and sd_likelihood, simulate 
%adaptive Bayes model (see Verstynen 2011) to generate predicted movement directions 
%theta_map. Also generate predicted time courses of subject's prior mean 
%and SD during experiment. Note that theta can contain multiple rows, each 
%being treated as an independent experiment. If 'sensory_uncertainty_flag'
%is 1, then the simulation is run using noisy target estimation. If it is
%0, the expected (noiseless) target location is used instead, which is
%necessary when fitting Bayesian parameters to the data.

%make sure prior learning rate 'B' and sd_likelihood are the right length 
%(equal to the number of trials to be simulated or scalars) 
if ~(isscalar(B) || (size(B,1) == size(theta,1) && size(B,2) == size(theta,2))) || ~(isscalar(sd_likelihood) || (size(sd_likelihood,1) == size(theta,1) && size(sd_likelihood,2) == size(theta,2)))
    error('B and sd_likelihood must be scalars (if constant) or the same size as theta (if changing over trials or simulations)')
end

if isscalar(B)
    B = repmat(B,size(theta,1),size(theta,2));
end

if isscalar(sd_likelihood)
    sd_likelihood = repmat(sd_likelihood,size(theta,1),size(theta,2));
end

%generate sensory signals x based on likelihood given true targets theta
x = zeros(size(theta));
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if sensory_uncertainty_flag == 1
            x(i,j) = theta(i,j) + sd_likelihood(i,j).*randn;
        else
            x(i,j) = theta(i,j);
        end
    end
end

%calculate prior updates and resulting MAP direction (i.e. the artificial 
%data) for each trial based on eqn 3 of Verstynen 2011
m_prior = zeros(size(theta,1),size(theta,2)+1);
m_prior(:,1) = m_prior_init;
sd_prior = zeros(size(theta,1),size(theta,2)+1);
sd_prior(:,1) = sd_prior_init;
theta_map = zeros(size(theta,1),size(theta,2));
for i = 1:size(x,1)
    for j = 1:size(x,2)
        sd_posterior = sqrt( ( (sd_prior(i,j).^-2) + (sd_likelihood(i,j).^-2) ).^-1 );
        if sd_prior(i,j) == 0
            theta_map(i,j) = m_prior(i,j);
        else
            %x must be expressed in range [m_prior-180,m_prior+180) to avoid circular effects
            x_adj = m_prior(i,j) + circDif(x(i,j),m_prior(i,j));
            theta_map(i,j) = ( (sd_posterior.^2)./(sd_prior(i,j).^2) ).*m_prior(i,j) + ( (sd_posterior.^2)./(sd_likelihood(i,j).^2) ).*x_adj;
            %make sure theta_map is in range [0,360)
            theta_map(i,j) = mod(theta_map(i,j),360);
        end
        %nonlinear function of distance between target and prior mean used for prior update
        d = circDif(theta(i,j),m_prior(i,j));
        trans_dist = d;
        %adjust prior mean based on circular dif between target and prior
        %mean to avoid circular effects
        m_prior(i,j+1) = m_prior(i,j) + B(i,j).*trans_dist;
        %make sure prior stays in [0,360)
        m_prior(i,j+1) = mod(m_prior(i,j+1),360); 
        sd_prior(i,j+1) = sqrt( (1-B(i,j)).*(sd_prior(i,j).^2) + B(i,j).*(trans_dist.^2) );
    end
end
