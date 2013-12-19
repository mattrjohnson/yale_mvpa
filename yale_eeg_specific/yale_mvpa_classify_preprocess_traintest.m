function [subs erps freqs] = mrj_eeg_classify_preprocess_traintest( paths, erps, freqs, start_time )

% very basic argument check before we get too far
if ~erps.use && ~freqs.use
    error('Either erps.use or freqs.use (or both) must be set to 1!');
end

n_subs =                                length( paths.sub_list );
n_trainsets =                           length( paths.train.p_sgc );
n_testsets =                            length( paths.test.p_sgc );

%---------------------------
% DO ACTUAL PRE-PROCESSING
%---------------------------

% pre-allocate struct
subs(n_subs).trainset(n_trainsets).erp_gav =                {};
subs(n_subs).trainset(n_trainsets).freq_gav =               {};
subs(n_subs).trainset(n_trainsets).n_conds =                [];
subs(n_subs).trainset(n_trainsets).n_trials_per_cond =      [];
subs(n_subs).trainset(n_trainsets).cond_names =             {};
subs(n_subs).trainset(n_trainsets).eegfile_inds =           {};
subs(n_subs).trainset(n_trainsets).trial_inds =             {};
subs(n_subs).trainset(n_trainsets).classifier_data =        {};

subs(n_subs).testset(n_testsets).erp_gav =                  {};
subs(n_subs).testset(n_testsets).freq_gav =                 {};
subs(n_subs).testset(n_testsets).n_conds =                  [];
subs(n_subs).testset(n_testsets).n_trials_per_cond =        [];
subs(n_subs).testset(n_testsets).cond_names =               {};
subs(n_subs).testset(n_testsets).eegfile_inds =             {};
subs(n_subs).testset(n_testsets).trial_inds =               {};
subs(n_subs).testset(n_testsets).classifier_data =          {};

for i=1:n_subs
    for j=1:n_trainsets
        % progress update
        disp(' ');
        disp(' ');
        disp(['Pre-processing subject ' int2str(i) ' of ' int2str(n_subs)]);
        disp([' - training set ' int2str(j) ' of ' int2str(n_trainsets)]);
        disp(['Script started: ' start_time.str]);
        disp(['Current time:   ' datestr(now)]);
        toc(start_time.tictoc);
        disp(' ');

        [subs(i).trainset(j) tp_inds] = mrj_eeg_classify_preprocess_oneset( paths, paths.train, erps, freqs, i, j );
        
        if erps.use && i==1 && j==1 %assume all subs/sets/conds/trials have same timepoints
            erps.tp_inds =              tp_inds;
        end
    end
    
    for j=1:n_testsets
        % progress update
        disp(' ');
        disp(' ');
        disp(['Pre-processing subject ' int2str(i) ' of ' int2str(n_subs)]);
        disp([' - test set ' int2str(j) ' of ' int2str(n_testsets)]);
        disp(['Script started: ' start_time.str]);
        disp(['Current time:   ' datestr(now)]);
        toc(start_time.tictoc);
        disp(' ');

        subs(i).testset(j) =           mrj_eeg_classify_preprocess_oneset( paths, paths.test, erps, freqs, i, j );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------
function [this_set tp_inds] = mrj_eeg_classify_preprocess_oneset( paths, paths_t, erps, freqs, subnum, setnum )

% reformat paths input to fool preprocessing functions
paths.eeg_subdir =                      paths_t.eeg_subdir;
paths.p_sgc =                           paths_t.p_sgc{ setnum };

% run MIPAVG3 to get per-trial ERPs
erp_data =                              mrj_eeg_classify_preprocess_domipavg( paths, erps, subnum );
this_set.erp_gav =                      mrj_eeg_classify_preprocess_erpgav( erp_data ); 
                                                            % calculate a grand average for posterity
    
% if using time-frequency info, do those calculations
if freqs.use
    freq_data =                         mrj_eeg_classify_preprocess_dofreqanal( erp_data, freqs );
    this_set.freq_gav =                 mrj_eeg_classify_preprocess_freqgav( freq_data );
end

% initialize some basic values; some don't really need to be redefined for every subject/set, but shouldn't hurt
this_set.n_conds =                      length(erp_data);
this_set.n_trials_per_cond =            zeros( this_set.n_conds, 1 );
for i=1:this_set.n_conds
    this_set.n_trials_per_cond(i) =     length(erp_data(i).trial);
    this_set.cond_names{i} =            erp_data(i).mrj_condition_name;
    this_set.eegfile_inds{i} =          erp_data(i).eegfile;
    this_set.trial_inds{i} =            erp_data(i).inds;
end

% update with number of trials
disp( 'Number of trials per condition: ' );
disp( this_set.n_trials_per_cond );

% prepare data for classification
tp_inds =                               mrj_get_tp_inds( erp_data(1).time{1}, erps.timepoints_to_use ); % so we can return to calling function
erps.tp_inds =                          tp_inds; % so we can pass to reshaping functions below

for i=1:this_set.n_conds
    if erps.use && freqs.use
        if erps.by_timepoint
            % first reshape individually
            erps_tmp =                  mrj_reshape_erp_for_classifier( erp_data, erps, i, this_set.n_trials_per_cond(i) );
            freqs_tmp =                 mrj_reshape_freq_for_classifier( freq_data, freqs, i );
            
            % check to make sure numbers of timepoints match
            if length(erps_tmp) ~= length(freqs_tmp)
                error('If using both ERPs and freqs, with by_timepoint == 1, the number of timepoints must be the same for both ERPs and freqs');
            end
            
            % now concatenate
            n_timepoints =              length( erps_tmp );
            erpfreq_tmp =               cell( n_timepoints, 1 );
            for j=1:n_timepoints
                erpfreq_tmp{j} =        [ erps_tmp{j}; freqs_tmp{j} ];
            end
            this_set.classifier_data{i} =    erpfreq_tmp;
        else
            this_set.classifier_data{i} =   [   mrj_reshape_erp_for_classifier( erp_data, erps, i, this_set.n_trials_per_cond(i) );
                                                mrj_reshape_freq_for_classifier( freq_data, freqs, i ) ];
        end
    elseif erps.use
        this_set.classifier_data{i} =    mrj_reshape_erp_for_classifier( erp_data, erps, i, this_set.n_trials_per_cond(i) );
    elseif freqs.use
        this_set.classifier_data{i} =    mrj_reshape_freq_for_classifier( freq_data, freqs, i );
    else
        error( 'Should never get here!' );
    end
end


%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------
function classifier_data = mrj_reshape_erp_for_classifier( erp_data, erps, cond_num, n_trials )

%make trial first dim, time second dim, remaining dims can be whatever depending on data type
classifier_data = zeros([n_trials fliplr(size(erp_data(cond_num).trial{1}))]);

for i=1:n_trials
    classifier_data(i,:,:) =            (erp_data(cond_num).trial{i})'; 
                                                            % classifier_data indexed by trial, timepoint, channel
end

if erps.convert_z %convert to z-scores
    for i=1:size(classifier_data,1) %loop thru trials
        for j=1:size(classifier_data,3) %loop thru channels
            classifier_data(i,:,j) =    mrj_zscores(classifier_data(i,:,j));
        end
    end
end
classifier_data =                       mrj_downsample( classifier_data(:,erps.tp_inds,erps.chans_to_use), erps.points_per_bin );
                                                            % classifier_data still indexed by trial, timepoint, channel
if erps.by_timepoint % make cell array by timepoint
    n_timepoints =                      size( classifier_data,2 );
    classifier_data_tmp =               cell( n_timepoints, 1 );
    for i=1:n_timepoints
        classifier_data_tmp{i} =        mrj_remove_nonfinite_features( squeeze( classifier_data(:,i,:) )' );
                                                            % each cell should end up n_channels (aka features) by n_trials
    end
    classifier_data =                   classifier_data_tmp;
else
    classifier_data =                   reshape( classifier_data, [size(classifier_data,1) numel(classifier_data)/size(classifier_data,1)])'; 
                                                                % now nfeatures by ntrials
    classifier_data =                   mrj_remove_nonfinite_features( classifier_data );
end


%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
function classifier_data = mrj_reshape_freq_for_classifier( freq_data, freqs, cond_num )

freq_data =                             freq_data{cond_num};
% freq_data starts out indexed by trial,channel,freq,timepoint
if freqs.by_timepoint % make cell array by timepoint
    n_trials =                          size( freq_data, 1 );
    n_timepoints =                      size( freq_data, 4 );
    n_features =                        size( freq_data, 2 ) * size( freq_data, 3 );
                                                            % multiply n_channels * n_freqs to get # of features at each timepoint
    classifier_data =                   cell( n_timepoints, 1 );
    for i=1:n_timepoints
        classifier_data{i} =            mrj_remove_nonfinite_features( reshape( squeeze( freq_data(:,:,:,i) ), [n_trials n_features] )' );
                                                            % each cell should end up n_features by n_trials
    end
else
    classifier_data =                   (reshape(freq_data,[size(freq_data,1) numel(freq_data)/size(freq_data,1)]))'; 
                                                            % now nfeatures by ntrials
    classifier_data =                   mrj_remove_nonfinite_features( classifier_data );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTILITY FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------
%---------------------------------------------
function tp_inds=mrj_get_tp_inds(tps,tp_range)

tp_inds=find(tps>=tp_range(1) & tps<=tp_range(2));


%-----------------------------------
%-----------------------------------
function outdata=mrj_zscores(indata)

outdata=(indata - mean(indata))/std(indata);


%-----------------------------------------------------
%-----------------------------------------------------
function outdata=mrj_downsample(indata,points_per_bin)
%assume first dimension is trials, second is time, can be as many other dimensions as you want

orig_dims=size(indata);
if length(orig_dims)>3
    indata=reshape(indata,[orig_dims(1:2) prod(orig_dims(3:end))]);
end %indata should now be 3 dimensions (or possibly 2)
    
while mod(size(indata,2),points_per_bin)
    indata=indata(:,1:end-1,:);
end
outdata_size=size(indata);
outdata_size(2)=outdata_size(2)/points_per_bin;
outdata=zeros(outdata_size);
for i=1:points_per_bin
    outdata=outdata+indata(:,i:points_per_bin:end,:);
end
outdata=outdata / points_per_bin;

if length(orig_dims)>3
    outdata=reshape(outdata,[size(outdata,1) size(outdata,2) orig_dims(3:end)]);
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function classifier_data = mrj_remove_nonfinite_features( classifier_data )

% classifier_data should be indexed nfeatures x ntrials at this point
nonfinite_features = any(~isfinite(classifier_data),2);
if any(nonfinite_features)
    classifier_data = classifier_data(~nonfinite_features,:);
end


%----------------------------------------------------------------
%----------------------------------------------------------------
function erp_gav = mrj_eeg_classify_preprocess_erpgav( erp_data )

% erp_gav output should be n_conds x 1 cell array, each containing an nchannels x ntimepoints matrix

n_conds =                               length( erp_data );
erp_gav =                               cell( n_conds, 1 );

for i=1:n_conds % loop over conditions
    n_trials =                          length( erp_data(i).trial );
    this_cond_sum =                     zeros( size(erp_data(i).trial{1}) );
    
    for j=1:n_trials % loop over individual trials
        this_cond_sum =                 this_cond_sum + erp_data(i).trial{j};
    end
    
    erp_gav{i} =                        this_cond_sum ./ n_trials;
end


%-------------------------------------------------------------------
%-------------------------------------------------------------------
function freq_gav = mrj_eeg_classify_preprocess_freqgav( freq_data )

% freq_gav output should be n_conds x 1 cell array, each containing an nchannels x nfrequencies x ntimepoints matrix

n_conds =                               length( freq_data );
freq_gav =                              cell( n_conds, 1 );

for i=1:n_conds % loop over conditions
    freq_size =                         size( freq_data{i} );
    n_trials =                          freq_size(1);
    this_cond_sum =                     zeros( freq_size(2:4) );
    
    for j=1:n_trials % loop over individual trials
        this_cond_sum =                 this_cond_sum + squeeze( freq_data{i}(j,:,:,:) );
    end
    
    freq_gav{i} =                       this_cond_sum ./ n_trials;
end
