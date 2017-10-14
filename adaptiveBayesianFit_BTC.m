function [sd_likelihood,B,ctxt_sd,subjects,fval,exitflag,subject_block_ct] = adaptiveBayesianFit_BTC(datapath,day,fb_type,block_list,m_prior_init,sd_prior_init,algo,varargin)
%For experiment data contained in folder datapath and with block names 
%specified by structure 'block_list', find likelihood SD and prior
%learning rate 'B' for adaptive Bayesian model as described in Verstynen 
%2011. 'day' is a string that specifies the days to be analyzed, and can be
%'Day1','Day2', or 'BothDays'; 'block_list' should have a field for each 
%day to be analyzed that specifies, in order, the blocks to be analyzed for 
%that day. 'algo' should be a structure specifying the fitting algorithm to 
%use in fmincon when fitting the model to each subject (i.e. fieldnames 
%should be subejcts, e.g. 'AB68', and values should be the algorithm name). 
%varargin should specify whether fitting is per-block (default) 
%or per-subject (all blocks for a subject fit at once). Output 'ctxt_sd' 
%specifies, for each subject, the context SDs presented to the subject on 
%each block (again, in the block order that blocks appear in 'block_list').
%'subject_block_ct' specifies the number of blocks that were analyzed for
%each subject (to allow verification that the correct number were
%analyzed). 'm_prior_init' and 'sd_prior_init' can be scalars if all
%subjects are to be fit based on the same initial prior mean and/or sd, or 
%can be vectors if each subject is to be fit with a different initial prior
%mean or sd (note this requires a priori knowledge of how many subjects'
%data will be fit).

if isempty(varargin) || strcmpi(varargin{1},'block')
%%%%NEEDS UPDATING
    error('By-block fitting not implemented')
%     %find list of appropriate subjects to analyze
%     [files,subjects] = findDataFiles(datapath,day,fb_type);
% 
%     current_dir = cd;
%     cd(datapath)
% 
%     %for each subject and block that matches a block specified in 
%     %'block_list', fit the model and record the context SD of the block
%     sd_likelihood = zeros(length(subjects),length(ctxt_sd_fit));
%     B = zeros(size(sd_likelihood));
%     fval = zeros(size(B));
%     exitflag = zeros(size(B));
%     ctxt_sd = zeros(size(B));
%     subject_block_flags = zeros(length(subjects),1); %keep track of number of appropriate blocks already loaded for each subject
%     for i = 1:length(files)
%         load(files{i},'exper');
%         if max(ctxt_sd_fit == exper.ctxt_sd) == 1
%             current_subject = getDataFileInfo(files{i});
%             row = find(subjects == current_subject);
%             subject_block_flags(row) = subject_block_flags(row) + 1;
%             col = subject_block_flags(row);
%             ctxt_sd(row,col) = exper.ctxt_sd;
%             load(files{i},'kdata');
%             [~,initang_data,targets] = initangError(kdata,0,0);
%             %fill in any missing data in initang_data by adding the mean error of
%             %the two surrounding trials to the target on the missing trial (or simply
%             %the error from the single preceding or following trial if the missing
%             %data is at the beginning or end of the block)
%             for j = 1:length(initang_data)
%                 if isnan(initang_data(j))
%                     id  = find(kdata(:,1) == j,1);
%                     targets(j) = cartesian2deg([kdata(id,7) kdata(id,8)]);
%                     tmp_err = nan(2,1);
%                     for k = (j-1):-1:1
%                         if ~isnan(initang_data(k))
%                             tmp_err(1) = circDif(initang_data(k),targets(k));
%                             break
%                         end
%                     end
%                     for k = (j+1):length(initang_data)
%                         if ~isnan(initang_data(k))
%                             tmp_err(2) = circDif(initang_data(k),targets(k));
%                             break
%                         end
%                     end
%                     tmp_err = mean_ignorenan(tmp_err,1);
%                     initang_data(j) = mod(targets(j) + tmp_err,360);
%                 end
%             end
%             clear tmp_err id
%             [sd_likelihood(row,col),B(row,col),fval(row,col),exitflag(row,col)] = fitBayesAdaptive(m_prior_init,sd_prior_init,targets,initang_data,algo.(['AB' num2str(current_subject)]));
%         end
%     end
% 
%     cd(current_dir)
    
elseif strcmpi(varargin{1},'subject')
    
    %find list of appropriate subjects to analyze
    if strcmpi(day,'Day1') || strcmpi(day,'Day2')
        [files,subjects] = findDataFiles(datapath,day,fb_type);
    elseif strcmpi(day,'BothDays')
        [files,subjects] = findDataFiles(datapath,'Day1',fb_type);
        files2 = findDataFiles(datapath,'Day2',fb_type);
        files = {files{:},files2{:}}; 
        clear files2
    end

    current_dir = cd;
    cd(datapath)
    
    %check if all subjects should share the same initial prior param values
    if length(m_prior_init) == 1
        m_prior_init = repmat(m_prior_init,1,length(subjects));
    end
    if length(sd_prior_init) == 1
        sd_prior_init = repmat(sd_prior_init,1,length(subjects));
    end

    %for each subject's block that matches a block specified in 'block_list' 
    %for the specified day, fit the model and record the context SD
    sd_likelihood = zeros(length(subjects),1);
    B = zeros(size(sd_likelihood));
    fval = zeros(size(B));
    exitflag = zeros(size(B));
    if strcmpi(day,'Day1')
        num_blocks = length(block_list.day1);
    elseif strcmpi(day,'Day2')
        num_blocks = length(block_list.day2);
    elseif strcmpi(day,'BothDays')
        num_blocks = length(block_list.day1) + length(block_list.day2);
    end
    default_block_length = 300; %this will be used to allocate space in the 'targets' and 'initang_data' array for data from each block. it must be equal to or larger than the number of data points contained in any single block for a single subject, otherwise data from some blocks will overlap data from other blocks in 'targets' and 'initang_data'
    ctxt_sd = zeros(length(subjects),num_blocks);
    subject_block_ct = zeros(length(subjects),1); %keep track of number of appropriate blocks already loaded for each subject
    targets = nan(length(subjects),num_blocks*default_block_length);
    initang_data = nan(length(subjects),num_blocks*default_block_length);
    for i = 1:length(files)
        [current_subject,block,day_number] = getDataFileInfo(files{i});
        col = cellfun(@(c)strcmpi(c,block),block_list.(lower(day_number))); %check if current file's block appears in the list of blocks to be analyzed, and in what place order
        if max(col) == 1
            load(files{i},'exper','kdata');
            block_length = length(exper.taskschedule);
            if block_length > default_block_length
                error('Default block length is too small')
            end
            row = find(subjects == current_subject);
            col = find(col == 1);
            if strcmpi(day_number,'Day2') && strcmpi(day,'BothDays')
                col = col + length(block_list.day1);
            end
            subject_block_ct(row) = subject_block_ct(row) + 1;
            ctxt_sd(row,col) = exper.ctxt_sd;
            [~,initang_data(row,(default_block_length*(col-1)+1):(default_block_length*(col-1)+block_length)),targets(row,(default_block_length*(col-1)+1):(default_block_length*(col-1)+block_length))]...
                = initangError(kdata,0,0);
            %fill in any missing data in initang_data by adding the mean error of
            %the two surrounding trials to the target on the missing trial (or simply
            %the error from the single preceding or following trial if the missing
            %data is at the beginning or end of the block)
            for j = 1:block_length
                if isnan(initang_data(row,default_block_length*(col-1)+j))
                    id  = find(kdata(:,1) == j,1);
                    targets(row,default_block_length*(col-1)+j) = cartesian2deg([kdata(id,7) kdata(id,8)]);
                    tmp_err = nan(2,1);
                    for k = (j-1):-1:1
                        if ~isnan(initang_data(row,default_block_length*(col-1)+k))
                            tmp_err(1) = circDif(initang_data(row,default_block_length*(col-1)+k),targets(row,default_block_length*(col-1)+k));
                            break
                        end
                    end
                    for k = (j+1):default_block_length
                        if ~isnan(initang_data(row,default_block_length*(col-1)+k))
                            tmp_err(2) = circDif(initang_data(row,default_block_length*(col-1)+k),targets(row,default_block_length*(col-1)+k));
                            break
                        end
                    end
                    tmp_err = mean_ignorenan(tmp_err,1);
                    initang_data(row,default_block_length*(col-1)+j) = mod(targets(row,default_block_length*(col-1)+j) + tmp_err,360);
                end
            end
            clear tmp_err id
        end
    end
    id = sum(isnan(targets)) == length(subjects);
    targets(:,id) = [];
    initang_data(:,id) = [];
    
    for i = 1:length(subjects)
        [sd_likelihood(i),B(i),fval(i),exitflag(i)] = fitBayesAdaptive(m_prior_init(i),sd_prior_init(i),targets(i,:),initang_data(i,:),algo.(['AB' num2str(current_subject)]));
    end

    cd(current_dir)
    
else
    
    error('Specification of by-block or by-subject fitting is unrecognized')
    
end
        

