function [simulated_VBias, mean_bias] = simulateVBias(subjects,m_prior_init,sd_prior_init,B_struct,sd_likelihood_struct,plot_flag)
%Based on the parameter values in 'm_prior_init', 'sd_prior_init',
%'B_struct', and 'sd_likelihood_struct', simluate the Bias part of the 
%Verstynen experiment based on the unique context/probe schedules found 
%among subject numbers in 'subjects'. 'm_prior_init' and 'sd_prior_init' are
%scalar values used for initiating the prior in all simulations. 'B_struct'
%and 'sd_likelihood_struct' are stuctures for which each field contains the
%results of fitting the Bayesian adaptive model to a specific experiment
%(e.g. Verstynen Bias, Bias Time Course, etc.) Fieldnames must match for 
%both structures. The different parameter values within each field are the 
%different fits to specific subjects. Iff 'plot_flag' is 1, plots will be
%generated.

%Set variable values
datapath = 'F:\CogGph\ExperimentReachRotation\Data\Accuracy only scoring';
block_list.bias = {'Block1','Block2','Block3','Block4'};
ctxt_m = 150;
unique_probe_dir = [30 60 90];
sensory_uncertainty_flag = 0;

%Find target schedules for each subject number in 'subjects', and retain
%for simulation.
[theta, probe_id, ctxt_sd] = targetList(datapath,subjects(1),block_list,'V');
theta = repmat(theta,length(subjects),1);
probe_id = repmat(probe_id,length(subjects),1);
ctxt_sd = repmat(ctxt_sd,length(subjects),1);
for i = 2:length(subjects)
    [theta(i,:), probe_id(i,:), ctxt_sd(i,:)] = targetList(datapath,subjects(i),block_list,'V');
end

%For each field in 'B_struct', create a corresponding field in
%'simulated_VBias'. Each field will contain another structure, with fields
%for the simulated data, probe biases, the corresponding targets, context 
%SD and probe trial numbers, and the biases, targets and probe trials 
%numbers for each of the probe distances. Also determine which 
%B/sd_likelihood parameters are valid (not NaN).
fn = fieldnames(B_struct);
valid_params = cell(1,length(fn));
simulated_VBias = struct();
for i = 1:length(fn)
    valid_params{i} = find(~isnan(B_struct.(fn{i})));
    %if there are only n subjects in the Verstynen Bias experiment but m
    %valid parameter sets in current B_struct field, where m > n, only take 
    %the first n parameter sets for simulating the n subjects 
    if length(valid_params{i}) > length(subjects)
        valid_params{i} = valid_params{i}(1:length(subjects));
    end
    simulated_VBias.(fn{i}) = struct(...
        'm_prior',zeros(length(valid_params{i}),size(theta,2)+1),...
        'sd_prior',zeros(length(valid_params{i}),size(theta,2)+1),...
        'theta',zeros(length(valid_params{i}),size(theta,2)),...
        'probes',zeros(length(valid_params{i}),size(probe_id,2)),...
        'theta_map',zeros(length(valid_params{i}),size(theta,2)),...
        'bias',zeros(length(valid_params{i}),size(probe_id,2)),...
        'ctxt_sd',zeros(length(valid_params{i}),size(theta,2)),...
        'theta30',nan(length(valid_params{i}),size(probe_id,2)),...
        'probes30',nan(length(valid_params{i}),size(probe_id,2)),...
        'bias30',nan(length(valid_params{i}),size(probe_id,2)),...
        'ctxt_sd30',zeros(length(valid_params{i}),size(probe_id,2)),...
        'theta60',nan(length(valid_params{i}),size(probe_id,2)),...
        'probes60',nan(length(valid_params{i}),size(probe_id,2)),...
        'bias60',nan(length(valid_params{i}),size(probe_id,2)),...
        'ctxt_sd60',zeros(length(valid_params{i}),size(probe_id,2)),...
        'theta90',nan(length(valid_params{i}),size(probe_id,2)),...
        'probes90',nan(length(valid_params{i}),size(probe_id,2)),...
        'bias90',nan(length(valid_params{i}),size(probe_id,2)),...
        'ctxt_sd90',zeros(length(valid_params{i}),size(probe_id,2))...
        );
end

%Simulate VBias data for each set of given B and sd_likelihood parameters.
%Also find biases at probe trials, and biases and probes for different 
%probe distances separated into different fields.
for i = 1:length(fn)
    for k = 1:length(valid_params{i})
        B = B_struct.(fn{i})(valid_params{i}(k));
        sd_likelihood = sd_likelihood_struct.(fn{i})(valid_params{i}(k));
        simulated_VBias.(fn{i}).theta(k,:) = theta(valid_params{i}(k),:);
        simulated_VBias.(fn{i}).probes(k,:) = probe_id(valid_params{i}(k),:);
        [simulated_VBias.(fn{i}).theta_map(k,:),...
            simulated_VBias.(fn{i}).m_prior(k,:),...
            simulated_VBias.(fn{i}).sd_prior(k,:)] = simulateAdaptiveBayesExp(theta(valid_params{i}(k),:),m_prior_init,sd_prior_init,B,sd_likelihood,sensory_uncertainty_flag);
        simulated_VBias.(fn{i}).bias(k,:) = calcBias(ctxt_m,theta(valid_params{i}(k),probe_id(valid_params{i}(k),:)),simulated_VBias.(fn{i}).theta_map(k,probe_id(valid_params{i}(k),:)));
        simulated_VBias.(fn{i}).ctxt_sd(k,:) = ctxt_sd(valid_params{i}(k),:);
        probe_dir = abs(circDif(theta(valid_params{i}(k),probe_id(valid_params{i}(k),:)),ctxt_m));
        for m = 1:length(unique_probe_dir)
            id = probe_dir == unique_probe_dir(m);
            simulated_VBias.(fn{i}).(['theta' num2str(unique_probe_dir(m))])(k,:) = [theta(valid_params{i}(k),probe_id(valid_params{i}(k),id)) nan(1,length(probe_id(valid_params{i}(k),:))-sum(id))];
            simulated_VBias.(fn{i}).(['probes' num2str(unique_probe_dir(m))])(k,:) = [probe_id(valid_params{i}(k),id)  nan(1,length(probe_id(valid_params{i}(k),:))-sum(id))];
            simulated_VBias.(fn{i}).(['bias' num2str(unique_probe_dir(m))])(k,:) = [calcBias(ctxt_m,theta(valid_params{i}(k),probe_id(valid_params{i}(k),id)),simulated_VBias.(fn{i}).theta_map(k,probe_id(valid_params{i}(k),id))) nan(1,length(probe_id(valid_params{i}(k),:))-sum(id))];
            simulated_VBias.(fn{i}).(['ctxt_sd' num2str(unique_probe_dir(m))])(k,:) = [ctxt_sd(valid_params{i}(k),probe_id(valid_params{i}(k),id)) nan(1,length(probe_id(valid_params{i}(k),:))-sum(id))];
        end
    end
    %trim off excess columns for biases at specific probe distances
    for m = 1:length(unique_probe_dir)
        id = sum(isnan(simulated_VBias.(fn{i}).(['theta' num2str(unique_probe_dir(m))]))) == size(simulated_VBias.(fn{i}).(['theta' num2str(unique_probe_dir(m))]),1);
        simulated_VBias.(fn{i}).(['theta' num2str(unique_probe_dir(m))])(:,id) = [];
        id = sum(isnan(simulated_VBias.(fn{i}).(['probes' num2str(unique_probe_dir(m))]))) == size(simulated_VBias.(fn{i}).(['probes' num2str(unique_probe_dir(m))]),1);
        simulated_VBias.(fn{i}).(['probes' num2str(unique_probe_dir(m))])(:,id) = [];
        id = sum(isnan(simulated_VBias.(fn{i}).(['bias' num2str(unique_probe_dir(m))]))) == size(simulated_VBias.(fn{i}).(['bias' num2str(unique_probe_dir(m))]),1);
        simulated_VBias.(fn{i}).(['bias' num2str(unique_probe_dir(m))])(:,id) = [];
        id = sum(isnan(simulated_VBias.(fn{i}).(['ctxt_sd' num2str(unique_probe_dir(m))]))) == size(simulated_VBias.(fn{i}).(['ctxt_sd' num2str(unique_probe_dir(m))]),1);
        simulated_VBias.(fn{i}).(['ctxt_sd' num2str(unique_probe_dir(m))])(:,id) = [];
    end
end     

%plot predicted bias for each simulated subject (mean +/- SE for all probes 
%at given distance for single simulation). Plot each simulated subject's 
%data in a different figure, with different subplots for different context 
%SDs showing the bias at each probe distance.
unique_ctxt_sd = unique(ctxt_sd);
for i = 1:length(fn)
    mean_bias.(fn{i}) = zeros(length(unique_ctxt_sd),length(unique_probe_dir),length(valid_params{i}));
    if plot_flag == 1
        figure
    end
    for k = 1:length(unique_ctxt_sd)
        for m = 1:length(unique_probe_dir)
            for j = 1:length(valid_params{i})
                id = simulated_VBias.(fn{i}).(['ctxt_sd' num2str(unique_probe_dir(m))])(j,:) == unique_ctxt_sd(k);
                mean_bias.(fn{i})(k,m,j) = mean_ignorenan(simulated_VBias.(fn{i}).(['bias' num2str(unique_probe_dir(m))])(j,id),2);
            end
        end
        if plot_flag == 1
            subplot(ceil(length(unique_ctxt_sd)/2),2,k)
            plot(unique_probe_dir,squeeze(mean_bias.(fn{i})(k,:,:)))
            xlabel('Probe Distance')
            ylabel('Inward bias (deg)')
            title({['Simulated Verstynen Bias Experiment, Ctxt SD ' num2str(unique_ctxt_sd(k))],['Using Bayesian Parameters Fit From ' fn{i} ' Experimental Data']})
            box('off')
            ylim([-1 45])
            set(gca,'XTick',unique_probe_dir)
        end
    end
end

if plot_flag == 1
    %plot predicted mean bias (mean +/- SE of mean bias of each simulation)
    %across different simulations with the same context SD, and simulated 
    %with parameters arising from fitting the same experiment type.
    colors = ['r','g','b','m'];
    h = zeros(1,length(unique_ctxt_sd));
    for i = 1:length(fn)
        figure
        hold on
        M = mean(mean_bias.(fn{i}),3);
        SE = std(mean_bias.(fn{i}),0,3)./sqrt(size(mean_bias.(fn{i}),3));
        for k = 1:length(unique_ctxt_sd)
            h(k) = errorbar(unique_probe_dir,M(k,:),SE(k,:),colors(k),'LineWidth',2);
        end
        ylim([0 30])
        set(gca,'XTick',unique_probe_dir)
        legend(h,'0 deg Ctxt SD','7.5 deg Ctxt SD','15 deg Ctxt SD','Uniform Ctxt')
        xlabel('Probe Distance')
        ylabel('Inward bias (deg)')
        title({'Mean +/- SE of Mean Bias, Simulated Verstynen Bias Experiment Across Blocks',['Using Bayesian Parameters Fit From ' fn{i} ' Experimental Data']})
    end
end
