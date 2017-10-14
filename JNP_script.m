%The following cells run various analyses and simulations needed for
%the JNP manuscript

%%
%Parameters

subject_struct.b = [68, 69, 73, 79, 81, 89, 93, 94, 95, 96];
subject_struct.n = [80, 82, 83, 85, 86, 87, 88, 90, 91, 97];
m_prior_init = 135;
sd_prior_init = 110;
datapath = 'F:\CogGph\ExperimentReachRotation\Data';

%% Fit adaptive Bayesian model to data from BTC Day 1
% This replicates Figure 3B

fb_type = 'b';
subjects.(fb_type) = subject_struct.(fb_type);
block_list = struct('day1', {{'Baseline', 'Block1', 'Block2', 'Block3', 'Block4', 'Block5', 'Block6'}});

%Use 'SQP' algorithm to fit adaptive Bayesian model to data.
SQP_algos = struct();
for i = 1:length(subjects.(fb_type))
    SQP_algos.(['AB', num2str(subjects.(fb_type)(i))]) = 'sqp'; 
end
[SQP_sd_likelihood,SQP_B,~,~,~,SQP_exitflag] = adaptiveBayesianFit_BTC(datapath,'Day1',fb_type,block_list,m_prior_init,...
    sd_prior_init,SQP_algos,'subject');

B_struct = struct('BTCDay1', zeros(1, length(subjects.(fb_type))));
sd_likelihood_struct = struct('BTCDay1', zeros(1, length(subjects.(fb_type))));
for i = 1:length(SQP_exitflag)
    B_struct.('BTCDay1')(i) = SQP_B(i);
    sd_likelihood_struct.('BTCDay1')(i) = SQP_sd_likelihood(i);
end

%Simulate Day 1 data (includes simulated bias and prior variance over time)
[simulated_BTC_Day1, aligned_bias_Day1, all_trials_arr_Day1, trial_blocks_arr_Day1] = simulateBTCDay1(subjects.(fb_type),m_prior_init,...
    sd_prior_init,B_struct,sd_likelihood_struct,block_list.day1);

%Simulate Verstynen experiment data (includes simulated bias and prior variance over time)
[simulated_VBias, mean_Vbias] = simulateVBias(subjects.(fb_type),m_prior_init,sd_prior_init,B_struct,sd_likelihood_struct,1);

clear SQP_algos SQP_sd_likelihood SQP_B SQP_fval i


%% Bootstrap simulations to check if targets in Retention block drive any bias increase
%Perform 'n_boot' bootstrap simulations of start of Day 2 assuming a naive
%prior, to allow comparison with actual retention data. Use sampling with
%replacement to choose individual subjects' fitted params to use in
%simulations.
n_boot = 5000;

fn = fieldnames(B_struct);
if length(fn) ~= 1 || length(fieldnames(sd_likelihood_struct)) ~= 1
    error('Unexpected number of fields in B_struct or sd_likelihood_struct');
elseif strcmpi(fn, fieldnames(sd_likelihood_struct)) ~= 1
    error('Unmatched fieldnames in B_struct and sd_likelihood_struct')
end

%Create input variables to 'simulateBTCDay2' for bootstrapped simulations
%and perform simulations
boot_samp = randi(10, 1, n_boot);
B_struct_bootstrap = struct(fn{1}, B_struct.(fn{1})(boot_samp));
sd_likelihood_struct_bootstrap = struct(fn{1}, sd_likelihood_struct.(fn{1})(boot_samp));
subjects_bootstrap = subjects.(fb_type)(boot_samp);
[~, bootstrap_aligned_bias_Day2] = simulateBTCDay2(subjects_bootstrap,m_prior_init,...
    sd_prior_init,B_struct_bootstrap,sd_likelihood_struct_bootstrap,{'Block1'});
sim_retention_bias_bootstrap = bootstrap_aligned_bias_Day2.BTCDay1(:, 1:10);

%Calculate within-subject means for all bootstrapped simulations and then 
%find 95% confidence intervals.
m = mean(sim_retention_bias_bootstrap, 2);
conf_int = sort(m);
conf_int = [conf_int(round(0.025*n_boot)), conf_int(round(0.975*n_boot))];

clear B_struct_bootstrap sd_likelihood_struct_bootstrap subjects_bootstrap fn i

