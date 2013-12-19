function erps = mrj_eeg_classify_config_erps
%-------------------------------------------


erps.use =                              1;
erps.highpass_filter =                  'detrend';          % 'detrend':                    for linear detrending only
                                                            % non-zero value (e.g., 0.1):   to specify a filter in Hz (positive values only)
                                                            % 0 or empty:                   no high-pass frequency filtering
erps.chans_to_use =                     3:33;
erps.timepoints_to_use =                [0 1.5];            % 0 point is what we would see in EEGAD; i.e., timepoints can be negative
erps.points_per_bin =                   10;                  % at 500Hz acquisition, e.g., 10 points per bin = 20ms per bin = resample to 50Hz
erps.convert_z =                        0;                  % convert each trial to Z-scores before classifying?
erps.regress_out_channels =             1:2;                % if this field is defined and non-empty, regress out the signal from the 
                                                            % specified channels (e.g., HEOG/VEOG). Do not include these channels in 
                                                            % the chans_to_use field