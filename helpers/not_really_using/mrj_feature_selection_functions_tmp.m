%-----------------------------------------------------------------------------
function outdata=mrj_do_feature_selection_alldata( indata, feature_selection )

disp(['Doing feature selection (on all data): ' feature_selection.type;])
outdata=rmfield(indata,'classify');
if strcmp(feature_selection.type, 'anova') % do ANOVA
    %check to make sure number of features matches across conditions
    n_trials_per_cond=zeros(length(indata),1);
    n_features = size(indata(1).classify, 1);
    for i=1:length(indata)
        if size(indata(i).classify, 1) ~= n_features
            error('All conditions do not appear to use the same number of features... cannot use anova feature selection');
        end
        
        n_trials_per_cond(i) = size(indata(i).classify, 2);
    end
    anova_ps = zeros(n_features,1);
    anova_grouping = [];
    for i=1:length(n_trials_per_cond)
        anova_grouping((end+1):(end+n_trials_per_cond(i)))=i;
    end
    
    anova_data=[indata.classify]; %should be nfeatures x ntrials (total, summed across all conditions)
    for i=1:n_features
        anova_ps(i) = anova1(anova_data(i,:), anova_grouping, 'off');
    end
    
    for i=1:length(indata)
        outdata(i).classify=indata(i).classify(anova_ps<feature_selection.thresh,:);
    end
else
    error(['Feature selection method ' feature_selection.type 'unknown or not yet implemented']);
end

%--------------------------------------------------------------------------
function [outdata inds]=mrj_do_feature_selection_traindata( indata, conds, feature_selection )
% indata (and outdata) are nfeatures x ntrials (total)
% conds is nconds x ntrials (total)

disp(['Doing feature selection (on training data): ' feature_selection.type;])
if strcmp(feature_selection.type, 'anova') % do ANOVA
    anova_ps = zeros(size(indata,1),1);
    anova_grouping = (conds') * ((1:size(conds,1))');
    
    for i=1:size(indata,1)
        anova_ps(i) = anova1(indata(i,:), anova_grouping, 'off');
    end
    
    inds=find(anova_ps<feature_selection.thresh);
    outdata=indata(inds,:);
else
    error(['Feature selection method ' feature_selection.type 'unknown or not yet implemented']);
end

%-----------------------------------------------------------------------------
function outdata=mrj_do_feature_selection_alldata_balanced( indata, feature_selection )

disp(['Doing feature selection (on all data -- balanced): ' feature_selection.type;])
outdata=rmfield(indata,'classify');
if strcmp(feature_selection.type, 'anova') % do ANOVA
    %check to make sure number of features matches across conditions
    n_conds=length(indata);
    n_trials_per_cond=zeros(n_conds,1);
    n_features = size(indata(1).classify, 1);
    for i=1:n_conds
        if size(indata(i).classify, 1) ~= n_features
            error('All conditions do not appear to use the same number of features... cannot use anova feature selection');
        end
        
        n_trials_per_cond(i) = size(indata(i).classify, 2);
    end
    
    anova_ps = zeros(feature_selection.nits,n_features);
    trials_used_per_condition = ceil(min(n_trials_per_cond)/2); %use approx 50% of trials, equal numbers for all conditions
    anova_data=zeros(trials_used_per_condition,n_conds);
    inds=zeros(n_conds,trials_used_per_condition);
    
    %loop through iterations
    for i=1:feature_selection.nits
        disp([' - feature selection iteration ' int2str(i)]);
        for j=1:n_conds
            inds_tmp=randperm(n_trials_per_cond(j));
            inds(j,:)=inds_tmp(1:trials_used_per_condition);
        end
        
        for j=1:n_features
            for k=1:n_conds
                anova_data(:,k)=indata(k).classify(j,inds(k,:));
            end
            anova_ps(i,j)=anova1(anova_data, [], 'off');
        end
    end
    
    anova_ps=mean(anova_ps);
    
    for i=1:n_conds
        outdata(i).classify=indata(i).classify(anova_ps<feature_selection.thresh,:);
    end
else
    error(['Feature selection method ' feature_selection.type 'unknown or not yet implemented']);
end

%------------------------------------------------------------------------------------------------------
function [outdata inds]=mrj_do_feature_selection_traindata_resample( indata, conds, feature_selection )
% indata (and outdata) are nfeatures x ntrials (total)
% conds is nconds x ntrials (total)

disp(['Doing feature selection (on training data -- resampling): ' feature_selection.type;])
if strcmp(feature_selection.type, 'anova') % do ANOVA
    n_features = size(indata, 1);
    n_conds=size(conds,1);
    n_trials_per_cond=size(conds,2) / n_conds; %unpredictable behavior if it doesn't divide evenly, but it should
    trials_used_per_condition = ceil(n_trials_per_cond/2); %use approx 50% of trials
    anova_ps = zeros(feature_selection.nits,n_features);
    condition_trial_inds=zeros(n_conds,n_trials_per_cond);
    anova_data=zeros(trials_used_per_condition,n_conds);
    inds=zeros(n_conds,trials_used_per_condition);
    
    for i=1:n_conds
        condition_trial_inds(i,:)=find(conds(i,:));
    end
    
    %loop through iterations
    for i=1:feature_selection.nits
        disp([' - feature selection iteration ' int2str(i)]);
        for j=1:n_conds
            inds_tmp=randperm(n_trials_per_cond);
            inds(j,:)=condition_trial_inds(j,inds_tmp(1:trials_used_per_condition));
        end
        
        parfor j=1:n_features
            anova_data=zeros(trials_used_per_condition,n_conds);
            for k=1:n_conds
                anova_data(:,k)=indata(j,inds(k,:));
            end
            anova_ps(i,j)=anova1(anova_data, [], 'off');
        end
    end
    
    anova_ps=mean(anova_ps);
        
    inds=find(anova_ps<feature_selection.thresh);
    outdata=indata(inds,:);
else
    error(['Feature selection method ' feature_selection.type 'unknown or not yet implemented']);
end
