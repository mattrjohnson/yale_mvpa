function [freq_data tp_inds] = mrj_eeg_classify_preprocess_dofreqanal( erp_data, freqs )

n_conds = length( erp_data );

% these assignments done so parfor will run without complaining
ft_cfg =                                freqs.ft_cfg;
normalize_to_baseline =                 freqs.normalize_to_baseline;
convert_z =                             freqs.convert_z;

% set timepoints and channels of interest in FieldTrip config
ft_cfg.toi =                            min(erp_data(1).time{1}) : ft_cfg.time_window_delta : max(erp_data(1).time{1});
ft_cfg.channel =                        erp_data(1).label(freqs.chans_to_use);

tp_inds =                               mrj_get_tp_inds(ft_cfg.toi,freqs.timepoints_to_use); 

parfor i=1:n_conds
    % run actual FieldTrip analysis
    fa_output =                         freqanalysis( ft_cfg, erp_data(i) );
                                                            % fa_output.powspec comes out ntrials x nchan x nfreq x ntimepoints
    % pull out power spectrum, flip it up
    powspec_alltrials =                 flipdim( fa_output.powspctrm, 3 );
                                                            % still indexed trial,channel,freq,timepoint, but with low frequencies at "bottom"
    
    % find negative timepoints / NaNs
    is_negtime =                        ( fa_output.time<=0 );
    powspec_isfinite =                  isfinite( squeeze( powspec_alltrials( 1,1,:,: ) ) ); 
                                                            % powspec_isfinite indexed by freq, timepoint
    % if normalizing to baseline, get baseline inds
    if normalize_to_baseline
        powspec_baseline_inds =         false(size(powspec_isfinite));
        for j=1:size(powspec_baseline_inds,1)               % loop thru frequencies to get baseline inds for each row
            baseline_inds =             (powspec_isfinite(j,:) & is_negtime);
            if ~any(baseline_inds)
                baseline_inds =         find(powspec_isfinite(j,:));
                if ~isempty(baseline_inds)
                    baseline_inds =     baseline_inds(1);
                end
            end
            powspec_baseline_inds(j,baseline_inds) =    true;
        end
    end
    
    if convert_z == 1                                       % z-score each frequency band for each trial
        powspec_alltrials =             mrj_zscore_over_dim( powspec_alltrials, 4 );
    elseif convert_z == 2                                   % z-score over trials
        powspec_alltrials =             mrj_zscore_over_dim( powspec_alltrials, 1 );
    end
    
    if normalize_to_baseline
        if convert_z == 1                                   % subtract out baseline average
            for j=1:size(powspec_alltrials,1)                   % loop thru trials
                for k=1:size(powspec_alltrials,2)               % loop thru channels
                    for m=1:size(powspec_alltrials,3)           % loop thru frequencies
                        if any(powspec_baseline_inds(m,:))
                            powspec_alltrials(j,k,m,:) =    powspec_alltrials(j,k,m,:) - mean(powspec_alltrials(j,k,m,powspec_baseline_inds(m,:)));
                        end
                    end
                end
            end
        else                                                % divide out baseline average
            for j=1:size(powspec_alltrials,1)                   % loop thru trials
                for k=1:size(powspec_alltrials,2)               % loop thru channels
                    for m=1:size(powspec_alltrials,3)           % loop thru frequencies
                        if any(powspec_baseline_inds(m,:))
                            powspec_alltrials(j,k,m,:) =    powspec_alltrials(j,k,m,:) / mean(powspec_alltrials(j,k,m,powspec_baseline_inds(m,:)));
                        end
                    end
                end
            end
        end
    end
    
%     if normalize_to_baseline && convert_z                   % convert to Z-scores and subtract out baseline period average
%         for j=1:size(powspec_alltrials,1)                   % loop thru trials
%             for k=1:size(powspec_alltrials,2)               % loop thru channels
%                 for m=1:size(powspec_alltrials,3)           % loop thru frequencies
%                     thisrow =                           squeeze(powspec_alltrials(j,k,m,:));
%                     powspec_alltrials(j,k,m,:) =        (thisrow - mean(thisrow(powspec_isfinite(m,:)))) / std(thisrow(powspec_isfinite(m,:)));
%                     powspec_alltrials(j,k,m,:) =        powspec_alltrials(j,k,m,:) - mean(powspec_alltrials(j,k,m,powspec_baseline_inds(m,:)));
%                 end
%             end
%         end
%     elseif normalize_to_baseline                            % don't convert to Z-scores; divide all timepoints by baseline period average
%         for j=1:size(powspec_alltrials,1)                   % loop thru trials
%             for k=1:size(powspec_alltrials,2)               % loop thru channels
%                 for m=1:size(powspec_alltrials,3)           % loop thru frequencies
%                     powspec_alltrials(j,k,m,:) =        powspec_alltrials(j,k,m,:) / mean(powspec_alltrials(j,k,m,powspec_baseline_inds(m,:)));
%                 end
%             end
%         end
%     elseif convert_z                                        % convert to Z-scores, but don't normalize to baseline at all
%         for j=1:size(powspec_alltrials,1)                   % loop thru trials
%             for k=1:size(powspec_alltrials,2)               % loop thru channels
%                 for m=1:size(powspec_alltrials,3)           % loop thru frequencies
%                     thisrow =                           squeeze(powspec_alltrials(j,k,m,:));
%                     powspec_alltrials(j,k,m,:) =        (thisrow - mean(thisrow(powspec_isfinite(m,:)))) / std(thisrow(powspec_isfinite(m,:)));
%                 end
%             end
%         end
%     end
                                                            %will be different from erps.tp_inds b/c this is post-"resampling" by the time-frequency analysis
    freq_data{i} =                                      powspec_alltrials(:,:,:,tp_inds);
    % cells of freq_data still indexed by trial,channel,freq,timepoint -- will reshape when it comes time to classify
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTILITY FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------
%---------------------------------------------
function tp_inds=mrj_get_tp_inds(tps,tp_range)
    
tp_inds=find(tps>=tp_range(1) & tps<=tp_range(2));


%----------------------------------------------------
%----------------------------------------------------
function data = mrj_zscore_over_dim( data, dim )
% convert a matrix of arbitrary dimensionality to zscores over a given dimension

dims =                                  ndims(data);

% we'll permute the data so the dimension in question is at the end
needs_permuting =                       1;
if dim == ndims(data)
    needs_permuting =                   0;
end
    
% we'll also convert data with 3+ dimensions to 2D data
needs_reshaping =                       0;
if dims>2
    needs_reshaping =                   1;
end

% permute data if necessary
if needs_permuting
    permdims =                          [setdiff(1:dims, dim), dim];
    data =                              permute( data, permdims );
end

% reshape data if necessary
if needs_reshaping
    orig_size =                         size( data ); % note: actually "original" size after possible permuting
    data =                              reshape( data, [prod(orig_size(1:end-1)) orig_size(end)] );
end

% should now have a 2-D array where we want to zscore over the second dimension
% start by establishing which, if any, rows have NaNs or Infs -- we'll handle those separately for speed
nonfinite_rows = any(~isfinite(data),2);
finite_rows = find(~nonfinite_rows)';   % make sure to include the ' so we don't screw up the for loops
nonfinite_rows = find(nonfinite_rows)'; % ditto

% zscore all rows that have only finite values
for i=finite_rows
    data(i,:) =                         (data(i,:) - mean(data(i,:))) / std(data(i,:));
end

% now handle rows with nonfinite values
for i=nonfinite_rows
    f_inds =                            isfinite( data(i,:) );
    if any(f_inds) % not much we can do if no values are finite
        f_vals =                        data(i,f_inds);
        data(i,f_inds) =                (f_vals - mean(f_vals)) / std(f_vals);
    end
end

% unreshape data if necessary
if needs_reshaping
    data =                              reshape( data, orig_size );
end

% unpermute data if necessary
if needs_permuting
    data =                              ipermute( data, permdims);
end
