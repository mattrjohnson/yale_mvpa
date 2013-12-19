function freqs = mrj_eeg_classify_config_freqs
%---------------------------------------------
% NOTE: ERPs will generally be processed even if they are not used for classification, 
%  so some of these settings will depend on the ones in mrj_eeg_classify_config_erps 
%  even if erps.use is 0 (e.g., freqs.chans_to_use can be identical to or a subset 
%  of erps.chans_to_use, but not a superset)


freqs.use =                             1;
freqs.chans_to_use =                    3:33;               % will generally be the same as erps.chans_to_use, but allow them to be separate just in case
freqs.timepoints_to_use =               [0 1.5];            % same notation as erps.timepoints_to_use above
freqs.convert_z =                       1;                  % convert signal to z-scores?
                                                            % 0: don't convert to z-scores
                                                            % 1: z-score each frequency band for each trial
                                                            % 2: z-score across trials for each freq/timepoint
freqs.normalize_to_baseline =           1;                  % express all values in each frequency band as ratio of baseline mean? 
                                                            % n.b.: if freqs.convert_z == 1, baseline is actually subtracted, not divieded, out
freqs.frequencies =                     2:2:128;            % in Hz
freqs.time_window_delta =               .05;                % how far to slide window each time, in seconds; probably good to make it a multiple of your sample size
freqs.method =                          'mtmconvol';        % choices so far: 'mtmconvol' or 'wltconvol'
                                                            % mtmconvol: if set, must also set freqs.taper and freqs.time_window
                                                            % wltconvol: if set, must also set freqs.width and freqs.gwidth
freqs.taper =                           'hanning';          % 'hanning' or possibly others (only hanning has been checked to work so far)
freqs.time_window =                     (2:(6-2)/(length(freqs.frequencies)-1):6) ./ freqs.frequencies;
                                                            % can be a scalar (same time window for all frequencies) or vector (same size as freqs.frequencies)
                                                            % the setting above goes from a 2 cycle time window at the lowest frequency to a 6 cycle time window at the highest
% freqs.width =                           6;                  % FieldTrip default is 7
% freqs.gwidth =                          3;                  % FieldTrip default is 3

% semi-user-editable if you know what you're doing, but most of the good stuff is up above
freqs.ft_cfg =                          [];
freqs.ft_cfg.output =                   'pow';
freqs.ft_cfg.channel =                  'xxx';
freqs.ft_cfg.method =                   freqs.method;
if strcmp(freqs.method,'mtmconvol')
    freqs.ft_cfg.taper =                freqs.taper;
    if isscalar(freqs.time_window)
        freqs.time_window =             ones(size(freqs.ft_cfg.foi))*freqs.time_window;
    end
    freqs.ft_cfg.t_ftimwin =            freqs.time_window(:);
elseif strcmp(freqs.method,'wltconvol')
    freqs.ft_cfg.width =                freqs.width;
    freqs.ft_cfg.gwidth =               freqs.gwidth;
end
freqs.ft_cfg.foi =                      freqs.frequencies;
freqs.ft_cfg.toi =                      [];                 % will set this once we read in FieldTrip .mat file
freqs.ft_cfg.time_window_delta =        freqs.time_window_delta; 
                                                            % only item in freqs.ft_cfg not intended for FieldTrip; used to calculate above (toi)
freqs.ft_cfg.keeptrials =               'yes';              % new for classification; keeps individual trials
