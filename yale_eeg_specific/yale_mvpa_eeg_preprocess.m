function [subs erps freqs] = mrj_eeg_classify_preprocess( paths, erps, freqs, start_time )

% very basic argument check before we get too far
if ~erps.use && ~freqs.use
    error('Either erps.use or freqs.use (or both) must be set to 1!');
end

n_subs = length(paths.sub_list);

%---------------------------
% DO ACTUAL PRE-PROCESSING
%---------------------------

% pre-allocate struct
subs(n_subs).n_conds =                  [];
subs(n_subs).cond_names =               {};
subs(n_subs).n_trials_per_cond =        [];
subs(n_subs).classifier_data =          {};
subs(n_subs).erp_gav =                  {};
subs(n_subs).freq_gav =                 {};
subs(n_subs).eegfile_inds =             {};
subs(n_subs).trial_inds =               {};

for i=1:n_subs
    % progress update
    disp(' ');
    disp(' ');
    disp(['Pre-processing subject ' int2str(i) ' of ' int2str(n_subs)]);
    disp(['Script started: ' start_time.str]);
    disp(['Current time:   ' datestr(now)]);
    toc(start_time.tictoc);
    disp(' ');
    
    % run MIPAVG3 to get per-trial ERPs
    erp_data =                          mrj_eeg_classify_preprocess_domipavg( paths, erps, i );
    if isfield(erps,'regress_out_channels') && ~isempty(erps.regress_out_channels)
        erp_data =                      mrj_eeg_classify_preprocess_regress_out_chans( erp_data, erps );
    end
    subs(i).erp_gav =                   mrj_eeg_classify_preprocess_erpgav( erp_data ); 
                                                            % calculate a grand average for posterity
    
    % if using time-frequency info, do those calculations
    if freqs.use
        freq_data =                     mrj_eeg_classify_preprocess_dofreqanal( erp_data, freqs );
        subs(i).freq_gav =              mrj_eeg_classify_preprocess_freqgav( freq_data );
                                                            % calculate a grand average for posterity
    end
    
    % initialize some basic values; some don't really need to be redefined for every subject, but shouldn't hurt
    subs(i).n_conds =                   length(erp_data);
    subs(i).n_trials_per_cond =         zeros( subs(i).n_conds, 1 );
    for j=1:subs(i).n_conds
        subs(i).n_trials_per_cond(j) =  length(erp_data(j).trial);
        subs(i).cond_names{j} =         erp_data(j).mrj_condition_name;
        subs(i).eegfile_inds{j} =       erp_data(j).eegfile;
        subs(i).trial_inds{j} =         erp_data(j).inds;
    end
        
    % update with number of trials per subject
    disp( ['Subject ' int2str(i) ' of ' int2str(n_subs) ', number of trials per condition: '] );
    disp( subs(i).n_trials_per_cond );
    
    % prepare data for classification
    if erps.use && i==1 %assume all subs/conds/trials have same timepoints
        erps.tp_inds = mrj_get_tp_inds( erp_data(1).time{1}, erps.timepoints_to_use );
    end
    
    for j=1:subs(i).n_conds
        if erps.use && freqs.use
            if erps.by_timepoint
                % first reshape individually
                erps_tmp =                  mrj_reshape_erp_for_classifier( erp_data, erps, j, subs(i).n_trials_per_cond(j) );
                freqs_tmp =                 mrj_reshape_freq_for_classifier( freq_data, freqs, j );
                
                % check to make sure numbers of timepoints match
                if length(erps_tmp) ~= length(freqs_tmp)
                    error('If using both ERPs and freqs, with by_timepoint == 1, the number of timepoints must be the same for both ERPs and freqs');
                end
                
                % now concatenate
                n_timepoints =              length( erps_tmp );
                erpfreq_tmp =               cell( n_timepoints, 1 );
                for k=1:n_timepoints
                    erpfreq_tmp{k} =        [ erps_tmp{k}; freqs_tmp{k} ];
                end
                subs(i).classifier_data{j} =    erpfreq_tmp;
            else
                subs(i).classifier_data{j} = [  mrj_reshape_erp_for_classifier( erp_data, erps, j, subs(i).n_trials_per_cond(j) );
                                                mrj_reshape_freq_for_classifier( freq_data, freqs, j ) ];
            end
        elseif erps.use
            subs(i).classifier_data{j} =    mrj_reshape_erp_for_classifier( erp_data, erps, j, subs(i).n_trials_per_cond(j) );
        elseif freqs.use
            subs(i).classifier_data{j} =    mrj_reshape_freq_for_classifier( freq_data, freqs, j );
        else
            error( 'Should never get here!' );
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
function erp_data = mrj_eeg_classify_preprocess_regress_out_chans( erp_data, erps )
% erp_data will have exactly the same structure as it does when passed in, but with specified channels regressed out
%  n.b.: the channels to-be-regressed-out will not change, as they would just become flat lines anyway if they were

disp(['Regressing out channels: ' mat2str(erps.regress_out_channels)]);

n_conds =                               length( erp_data );

for i=1:n_conds % loop over conditions
    for j=1:length( erp_data(i).trial ) % loop over inidividual trials within condition
        this_trial_data =               erp_data(i).trial{j}; % indexed channel x timepoint
        x =                             this_trial_data( erps.regress_out_channels, : )'; % x is now timepoint by channel
        for k=erps.chans_to_use % loop over channels of interest
            [b bint r] =                regress( this_trial_data(k,:)', x ); %#ok<ASGLU>
            this_trial_data(k,:) =      r;
        end
        erp_data(i).trial{j} =          this_trial_data;
    end
end

disp(' ... done');