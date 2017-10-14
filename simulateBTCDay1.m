function [simulated_BTC, aligned_bias, all_trials_arr, trial_blocks_arr] = simulateBTCDay1(subjects,m_prior_init,sd_prior_init,B_struct,sd_likelihood_struct,varargin)
%Based on the parameter values in 'm_prior_init', 'sd_prior_init',
%'B_struct', and 'sd_likelihood_struct', simluate Day 1 of the Bias Time
%Course experiment based on the unique context/probe schedules found among
%subject numbers in 'subjects'. 'B_struct' and 'sd_likelihood_struct' are 
%structures for which each field contains the results of fitting the 
%Bayesian adaptive model to a specific experiment (e.g. Verstynen Bias, 
%Bias Time Course Day 1, etc.) Fieldnames must match for both structures. 
%The different parameter values within each field are the different fits to 
%specific subjects. 'm_prior_init' and 'sd_prior_init' are either scalar 
%values used for initiating the prior in all simulations, or are structures 
%with fields containing vectors of initial prior parameters (for each 
%subject) that are to be used when running simulations based on the 
%corresponding field/vector elements in B_struct and sd_likelihood_struct.
%'varargin' can contain a cell array of block names to be simulated. If not
%supplied, it is created in code.
%
%Output 'aligned_bias' contains all subjects' mean biases in each bin of
%probe trials, with bins the same for all subjects and created based on the
%most common probe schedule (same as during actual data analysis of mean 
%bias over time; see plotMeanBTCData). 'all_trials_arr' contains the
%corresponding initial trial number of each bin for easy plotting, while
%'trial_block_arr' contains the corresponding block name for each bin.

%Set variable values
if strncmpi(computer,'mac',3)
    dr = '/Volumes/Untitled';
elseif strncmpi(computer,'pc',2)
    dr = 'F:';
end
datapath = [dr '/CogGph/ExperimentReachRotation/Data'];
if isempty(varargin)
    block_list.day1 = {'Baseline','Block1','Block2','Block3','Block4','Block5','Block6'};
else
    block_list.day1 = varargin{1};
end
ctxt_m = 135;
sensory_uncertainty_flag = 0;

%If m_prior_init and sd_prior_init are scalars, transform into structures
%for use during simulations
fn = fieldnames(B_struct);
if ~isstruct(m_prior_init)
    tmp = m_prior_init;
    clear m_prior_init %to suppress warning about overwriting variable with a structure
    for i = 1:length(fn)
        m_prior_init.(fn{i}) = repmat(tmp,length(B_struct.(fn{i})),1);
    end
    clear tmp
end
if ~isstruct(sd_prior_init)
    tmp = sd_prior_init;
    clear sd_prior_init %to suppress warning about overwriting variable with a structure
    for i = 1:length(fn)
        sd_prior_init.(fn{i}) = repmat(tmp,length(B_struct.(fn{i})),1);
    end
    clear tmp
end

%Find target schedules for each subject number in 'subjects', and retain
%for simulation.
[~, ~, blocks] = targetList(datapath,subjects(1),block_list,'BTC');
theta = cell(length(subjects),1);
probe_id = cell(length(subjects),1);
blocks = cell(length(subjects),length(blocks));
for i = 1:length(subjects)
    [theta{i}, probe_id{i}, blocks(i,:)] = targetList(datapath,subjects(i),block_list,'BTC');
end

%For each field in 'B_struct', create a corresponding field in
%'simulated_BTC'. Each field will contain another structure, with fields
%for the simulated data, probe biases, and the corresponding targets and 
%probe trial numbers. Also determine which B/sd_likelihood parameters are 
%valid (not NaN).
valid_params = cell(1,length(fn));
simulated_BTC = struct();
for i = 1:length(fn)
    valid_params{i} = find(~isnan(B_struct.(fn{i})));
    %if there are only n subjects in the Bias Time Course experiment but m
    %valid parameter sets in current B_struct field, where m > n, only take 
    %the first n parameter sets for simulating the n subjects 
    if length(valid_params{i}) > length(subjects)
        valid_params{i} = valid_params{i}(1:length(subjects));
    end
    simulated_BTC.(fn{i}) = struct(...
        'm_prior',{cell(length(valid_params{i}),1)},...
        'sd_prior',{cell(length(valid_params{i}),1)},...
        'theta',{cell(length(valid_params{i}),1)},...
        'probes',{cell(length(valid_params{i}),1)},...
        'theta_map',{cell(length(valid_params{i}),1)},...
        'bias',{cell(length(valid_params{i}),1)},...
        'blocks',{cell(length(valid_params{i}),1)},...
        'probe_blocks',{cell(length(valid_params{i}),1)}...
        );
end

%Simulate BTC data for each set of given B and sd_likelihood parameters.
%Also find biases at probe trials.
for i = 1:length(fn)
    for k = 1:length(valid_params{i})
        B = B_struct.(fn{i})(valid_params{i}(k));
        sd_likelihood = sd_likelihood_struct.(fn{i})(valid_params{i}(k));
        mpi = m_prior_init.(fn{i})(valid_params{i}(k));
        spi = sd_prior_init.(fn{i})(valid_params{i}(k));
        simulated_BTC.(fn{i}).theta{k} = theta{valid_params{i}(k)};
        simulated_BTC.(fn{i}).probes{k} = probe_id{valid_params{i}(k)};
        [simulated_BTC.(fn{i}).theta_map{k},...
            simulated_BTC.(fn{i}).m_prior{k},...
            simulated_BTC.(fn{i}).sd_prior{k}] = simulateAdaptiveBayesExp(theta{valid_params{i}(k)},mpi,spi,B,sd_likelihood,sensory_uncertainty_flag);
        simulated_BTC.(fn{i}).bias{k} = calcBias(ctxt_m,theta{valid_params{i}(k)}(probe_id{valid_params{i}(k)}(~isnan(probe_id{valid_params{i}(k)}))),...
            simulated_BTC.(fn{i}).theta_map{k}(probe_id{valid_params{i}(k)}(~isnan(probe_id{valid_params{i}(k)}))));
        simulated_BTC.(fn{i}).blocks{k} = blocks(valid_params{i}(k),:);
        simulated_BTC.(fn{i}).probe_blocks{k} = blocks(valid_params{i}(k),probe_id{valid_params{i}(k)}(~isnan(probe_id{valid_params{i}(k)})));
    end
end     

%plot predicted bias for each simulation
colors = ['r','g','b','m','c','y'];
for i = 1:length(fn)
    figure
    hold on
    for k = 1:length(valid_params{i})
        plot(probe_id{valid_params{i}(k)}(~isnan(probe_id{valid_params{i}(k)})),simulated_BTC.(fn{i}).bias{k},colors(rem(k,length(colors))+1))
    end
    xlabel('Day 1 Trial Number')
    ylabel('Inward bias (deg)')
    title({'Simulated Day 1 Bias Time Course Data Using',['Bayesian Parameters Fit From ' fn{i} ' Experimental Data']})
    box('off')
    ylim([-1 15])
end

%check for schedule similarity between subjects, create trial bin 
%boundaries for probe trials, and find avg bias for each subject on probe 
%trials within each bin for each subject, as well as record probe trial 
%positions for later plotting
num_probe_dif_thresh = 2;
for i = 1:length(fn)
    aligned_bias.(fn{i}) = cell(1,length(block_list.day1));
    all_trials_arr.(fn{i}) = cell(1,length(block_list.day1));
    trial_blocks_arr.(fn{i}) = {};
    for m = 1:length(block_list.day1)
        probe_check = zeros(length(valid_params{i}),20); 
        for j = 1:length(valid_params{i})
            tmp = probe_id{valid_params{i}(j)}(strcmpi(blocks(valid_params{i}(j),probe_id{valid_params{i}(j)}),block_list.day1{m}));
            probe_check(j,1:length(tmp)) = tmp; 
        end
        %find most common probe schedule for current block
        [ur,~,idu] = unique(probe_check,'rows');
        standard_sched = ur(mode(idu),:);
        standard_sched = standard_sched(standard_sched > 0);
        %check all other schedules do not deviate in number of probes
        %from standard schedule by more than 'num_probe_dif_thresh'
        num_standard_probes = sum(standard_sched > 0);
        for k = 1:size(ur,1)
            if abs(num_standard_probes - sum(ur(k,:) > 0)) > num_probe_dif_thresh
                error(['Large difference in number of probes for different subjects, Block ' block_list.day1{m}])
            end
        end
        %for each trial in standard probe schedule, place the name of the
        %current block in the overall list of blocks corresponding to each
        %standard schedule probe trial
        trial_blocks_arr.(fn{i})(length(trial_blocks_arr.(fn{i}))+1:length(trial_blocks_arr.(fn{i}))+num_standard_probes) = block_list.day1(m);
        %create probe trial bin boundaries and record probe trials for later
        %plotting
        trial_boundaries = [standard_sched(1) (standard_sched(1:end-1)+standard_sched(2:end))/2 standard_sched(end)];
        all_trials_arr.(fn{i}){m} = standard_sched;
        %for each bin defined by 'trial_boundaries', find probe trials within
        %the bin for each valid subject and create corresponding data point in
        %'aligned_bias' that is the mean of all probe trials in the bin for the
        %given subject
        aligned_bias.(fn{i}){m} = nan(length(valid_params{i}),length(trial_boundaries)-1);
        for j = 1:length(valid_params{i})
            tmp_probe_trials = probe_id{valid_params{i}(j)}(strcmpi(simulated_BTC.(fn{i}).probe_blocks{j},block_list.day1{m}));
            tmp_bias = simulated_BTC.(fn{i}).bias{j}(strcmpi(simulated_BTC.(fn{i}).probe_blocks{j},block_list.day1{m}));
            for k = 1:length(trial_boundaries)-1
                if k == length(trial_boundaries)-1
                    id = find(tmp_probe_trials >= trial_boundaries(k) & tmp_probe_trials <= trial_boundaries(k+1));
                else
                    id = find(tmp_probe_trials >= trial_boundaries(k) & tmp_probe_trials < trial_boundaries(k+1));
                end
                aligned_bias.(fn{i}){m}(j,k) = mean_ignorenan(tmp_bias(id),2);
            end
        end
        clear tmp ur idu standard_sched probe_check
    end
    aligned_bias.(fn{i}) = cell2mat(aligned_bias.(fn{i}));
    all_trials_arr.(fn{i}) = cell2mat(all_trials_arr.(fn{i}));
end

%plot predicted bias mean +/- SE across different simulations simulated 
%with  parameters arising from fitting the same experiment type.
for i = 1:length(fn)
    figure
    M = mean_ignorenan(aligned_bias.(fn{i}),1);
    SE = std_ignorenan(aligned_bias.(fn{i}),1)./sqrt(sum(~isnan(aligned_bias.(fn{i}))));
    hold on
    fillErrorbar(all_trials_arr.(fn{i}),M,SE,'b',[0.7 0.7 1]);
    ylim([0 8])
    xlabel('Day 1 Trial Number')
    ylabel('Inward bias (deg)')
    title({'Mean +/- SE Simulated Day 1 Bias Time Course Data Using',['Bayesian Parameters Fit From ' fn{i} ' Experimental Data']})
end

