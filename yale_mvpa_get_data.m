function [yale_mvpa_config, subs] = yale_mvpa_get_data( yale_mvpa_config )

% get data, either via pre-processing or loading from file
if yale_mvpa_config.files.load_preprocessed_data
    if yale_mvpa_config.files.save_preprocessed_data
        error('Cannot set both files.load_preprocessed_data and files.save_preprocessed_data to 1');
    end
    
    if yale_mvpa_config.files.load_data_via_script
        error('Cannot set both files.load_preprocessed_data and files.load_data_via_script to 1');
    end
    
    loaded_data = load( yale_mvpa_config.files.load_preprocessed_data_fname );
    if ~isfield(loaded_data,'subs')
        error(['Load from ' yale_mvpa_config.files.load_preprocessed_data_fname ' failed; ''subs'' variable not found.']);
    end
    subs = loaded_data.subs;
    
    if yale_mvpa_config.files.save_results
        if ~isfield(loaded_data,'erps') || isempty(loaded_data.erps)
            loaded_data.erps = 'Preprocessed data loaded from file. No ''erps'' variable found (or it was empty).';
        end

        if ~isfield(loaded_data,'freqs') || isempty(loaded_data.freqs)
            loaded_data.freqs = 'Preprocessed data loaded from file. No ''freqs'' variable found (or it was empty).';
        end
        
        disp(['Preprocessed data loaded from file. Note that ''yale_mvpa_config.erps'' and ''yale_mvpa_config.freqs'' in the ' ...
              'saved results file will reflect values loaded from this file, not whatever is specified in your configuration ' ...
              'scripts (if those happen to differ). If you are not using EEG data (i.e. modality is fMRI) or not doing a '     ...
              'time-frequency analysis as part of your preprocessing, then those settings may contain placeholder values.']);
        
        yale_mvpa_config.erps  = loaded_data.erps;
        yale_mvpa_config.freqs = loaded_data.freqs;
    end
    clear loaded_data; %save some memory, maybe
elseif yale_mvpa_config.files.load_data_via_script
    try
        subs = eval(yale_mvpa_config.files.load_data_via_script_mfile);
    catch %#ok<CTCH>
        error('Loading data via your specified script failed. Either yale_mvpa_config.files.load_data_via_script_mfile is empty, or the file could not be found, or it did not return any output.');
    end
elseif strcmp( yale_mvpa_config.general.modality, 'eeg' ) %deprecate later and fold into yale_mvpa_config.files.load_data_via_script
    [subs, yale_mvpa_config] = yale_mvpa_eeg_preprocess( yale_mvpa_config );
else
    error('No method of fMRI preprocessing implemented yet! Either change your modality setting or load data from a file for now.');
end

% possibly save pre-processed data
if yale_mvpa_config.files.save_preprocessed_data
    do_save = yale_mvpa_check_save_ok( yale_mvpa_config.files.save_preprocessed_data_fname );
    if do_save
        if strcmp( yale_mvpa_config.general.modality, 'eeg' )
            erps = yale_mvpa_config.erps; %#ok<NASGU>
            freqs = yale_mvpa_config.freqs; %#ok<NASGU>
            save( yale_mvpa_config.files.save_preprocessed_data_fname, 'subs', 'erps', 'freqs' );
            clear erps freqs;
        else
            error('No method of fMRI preprocessing, I said! You shouldn''t even be able to get to this line of code!');
        end
    end
end


