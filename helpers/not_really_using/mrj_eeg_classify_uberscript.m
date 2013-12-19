function mrj_eeg_classify_uberscript

%overall/mipavg settings
paths_function='mrj_eeg_paths_mac';
mipavg_files_folder='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files';

%for testing algorithms to ensure proper function
classifier.shuffle_data_randomly =0; %should be 0 most of the time; other values should yield chance performance
% n.b.: 1 shuffles features/trials randomly within condition (should maybe shuffle between conditions too to get truly uninformative data...)
%       2 pre-analyzes data as usual, but when time to classify comes, just sticks in random data of the same size

%general classification settings
classifier.kfold=5; %i.e. divide trials into fifths or whatnot for training/testing
classifier.nits=10; %number of times to iterate classifier

% classifier.trainfunc=@train_matlabsvm_mrj;
% classifier.testfunc=@test_matlabsvm_mrj;
% classifier.args=[]; %classifier.args will depend on what your train and test function are

classifier.trainfunc=@train_logreg;
classifier.testfunc=@test_logreg;
classifier.args.penalty=100; %classifier.args will depend on what your train and test function are

% classifier.trainfunc=@train_smlr;
% classifier.testfunc=@test_smlr;
% classifier.args.verbose=0; %classifier.args will depend on what your train and test function are

% classifier.trainfunc=@train_ridge;
% classifier.testfunc=@test_ridge;
% classifier.args.penalty=100; %classifier.args will depend on what your train and test function are

%settings for feature selection
feature_selection.use = 1;
feature_selection.type = 'anova';
feature_selection.thresh = .05;

%settings for classifying based on ERPs
erps.use=1;
erps.chans_to_use=3:33;
erps.timepoints_to_use=[0 1.5]; %0 point is what we would see in EEGAD; i.e., timepoints can be negative
erps.points_per_bin=10; %at 500Hz acquisition, e.g., 10 points per bin = 20ms per bin = resample to 50Hz
erps.convert_z=0; %convert each trial to Z-scores before classifying?

%settings for classifying based on time-frequency data
freqs.use=1;
freqs.chans_to_use=3:33; %will generally be the same as erps.chans_to_use, but allow them to be separate just in case
freqs.timepoints_to_use=[0 1.5]; %same notation as erps.timepoints_to_use above
freqs.convert_z=1; %convert signal in each frequency band to z-scores?
freqs.normalize_to_baseline=1; %express all values in each frequency band as ratio of baseline mean? (if both set, occurs after convert_to_zscore & baseline is actually subtracted, not divieded, out)
freqs.method='mtmconvol'; %at the moment this is the only option that will work
freqs.taper='hanning'; %ditto
freqs.frequencies=2:2:120; %in Hz
freqs.time_window_delta=.05; %how far to slide window each time, in seconds; probably good to make it a multiple of your sample size
freqs.time_window=(2:(6-2)/(length(freqs.frequencies)-1):6) ./ freqs.frequencies; %can be a scalar (same time window for all frequencies) or vector (same size as freqs.frequencies)
  % the setting above goes from a 2 cycle time window at the lowest frequency to a 6 cycle time window at the highest

%semi-user-editable if you know what you're doing, but most of the good stuff is up above
freqs.ft_cfg                  = [];
freqs.ft_cfg.output           = 'pow';
freqs.ft_cfg.channel          = 'xxx';
freqs.ft_cfg.method           = freqs.method;
freqs.ft_cfg.taper            = freqs.taper;
freqs.ft_cfg.foi              = freqs.frequencies;
if isscalar(freqs.time_window)
    freqs.time_window       = ones(size(freqs.ft_cfg.foi))*freqs.time_window;
end
freqs.ft_cfg.t_ftimwin        = freqs.time_window(:);
freqs.ft_cfg.toi              = []; %will set this once we read in FieldTrip .mat file
freqs.ft_cfg.time_window_delta= freqs.time_window_delta; %only item in freqs.ft_cfg not intended for FieldTrip; used to calculate above (toi)
freqs.ft_cfg.keeptrials       = 'yes'; %new for classification; keeps individual trials


%override GUI with these if desired

%for localizer
% sub_dirs=   [   '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100001/'
%                 '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100003/'
%                 '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100004/'
%                 '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100006/'
%                 '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100009/'
%             ];
% firstsub_eeg_dir='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100001/loc/';
% p_arf='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/refnegprime.arf';
% p_sgc='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/refnegprime_loc.sgc';
% p_mon='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/Cap34_ChinRef.mon';

%for refresh category
sub_dirs=   [   '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100001/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100002/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100003/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100004/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100006/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100007/'
                '/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100009/'
            ];
firstsub_eeg_dir='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/eeg_refnegprime4/20100001/maintask/';
p_arf='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/refnegprime.arf';
p_sgc='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/refnegprime_refresh_by_category.sgc';
p_mon='/Users/matt/Desktop/docs_tmp/work/refnegprime_tmp/standard_files/Cap34_ChinRef.mon';

%-------no user-editable parameters after this point

owd=pwd;


%===========%
% set paths %
%===========%
try
    eval(paths_function);
catch %#ok<CTCH>
    error('Could not evaluate paths function!');
end


%=============================================%
% get list of folders with .EEG files in them %
%=============================================%
if ~exist('sub_dirs','var') || isempty(sub_dirs)
    sub_dirs=spm_select(Inf,'dir','Select subject folders');
end
n_subs=size(sub_dirs,1);
cd(deblank(sub_dirs(1,:)));
if ~exist('firstsub_eeg_dir','var') || isempty(firstsub_eeg_dir)
    firstsub_eeg_dir=spm_select(1,'dir','Select folder with EEG files');
end
for i=1:n_subs %loop thru subjects
    eeg_dirs{i}=strrep(firstsub_eeg_dir,deblank(sub_dirs(1,:)),deblank(sub_dirs(i,:))); %#ok<AGROW>
end


%=======================================%
% get .arf, .sgc, .mon files for mipavg %
%=======================================%
if ~isempty(mipavg_files_folder)
    cd(mipavg_files_folder);
end
if ~exist('p_arf','var') || isempty(p_arf)
    p_arf=spm_select(1,  '\.arf','Select .arf file' );
end
if ~exist('p_sgc','var') || isempty(p_sgc)
    p_sgc=spm_select(1,  '\.sgc','Select .sgc file' );
end
if ~exist('p_mon','var') || isempty(p_mon)
    p_mon=spm_select(1,  '\.mon','Select .mon file' );
end


%=================================%
% high-pass filter or detrending? %
%=================================%
flags={'-triponly'}; %assume we don't care about saving averaged ERPs
btn=questdlg('Apply high-pass filter?',' ','Yes','Detrend only','No','No');
if strcmp(btn,'Yes')
    filt=spm_input('Cutoff frequency: ','','s');
    flags{end+1}=['-f' filt ',0'];
elseif strcmp(btn,'Detrend only')
    flags{end+1}='-lin';
end


%=================================
% keep track of how long it's been
%=================================
start_time=datestr(now);
tstart = tic;


%==========================================================================%
% loop through subjects, do full analysis on each before going to next sub %
%==========================================================================%
n_features_per_sub=zeros(n_subs,1);
for i=1:n_subs
    %================
    % progress update
    %================
    disp(' ');
    disp(' ');
    disp(['Analyzing subject ' int2str(i) ' of ' int2str(n_subs)]);
    disp(['Script started: ' start_time]);
    disp(['Current time:   ' datestr(now)]);
    toc(tstart);
    disp(' ');
    
    %====================================================================%
    % get .EEG files to work on (just skip this sub if they don't exist) %
    %====================================================================%
    try
        cd(deblank(eeg_dirs{i}));
    catch %#ok<CTCH>
        disp(['No such folder as ' deblank(eeg_dirs{i}) '; skipping']);
        continue;
    end
    eeg_files=cellstr(spm_select('List',pwd,'.EEG$'));
    
    %==================================================================%
    % generate a suitable temporary filename for fieldtrip-style input %
    %==================================================================%
    mip_output_tempname=mrj_tempname;
    while exist([mip_output_tempname '.avg'],'file') || exist([mip_output_tempname '.log'],'file') || exist([mip_output_tempname '.hdr'],'file') || exist([mip_output_tempname '_fieldtrip.mat'],'file')
        mip_output_tempname=mrj_tempname;
    end

    %========================================================%
    % run mipavg to pre-process trials and split into epochs %
    %========================================================%
    mrj_mipavg3(eeg_files{:},p_arf,p_sgc,p_mon,[mip_output_tempname '.avg'],[mip_output_tempname '.log'],flags{:});
    
    %=================================================%
    % load in trial epoch data, delete temporary file %
    %=================================================%
    erp_trialdata=load([mip_output_tempname '_fieldtrip.mat']);
    erp_trialdata=erp_trialdata.mipavg_fieldtrip_data;
    delete([mip_output_tempname '_fieldtrip.mat']);
    
    %========================%
    % define some key values %
    %========================%
    clear data; %make sure we start each subject with a clean slate
    n_conds=length(erp_trialdata); %some of these don't really need to be redefined for every subject, but shouldn't hurt
    n_trials_per_cond=zeros(n_conds,1);
    for j=1:n_conds
        n_trials_per_cond(j)=length(erp_trialdata(j).trial);
        data(j).cond_name=erp_trialdata(j).mrj_condition_name; %#ok<AGROW>
    end
    if i==1
        all_subs_trials_per_cond=zeros(n_subs,n_conds);
    end
    disp(['Subject ' int2str(i) ' of ' int2str(n_subs) ', number of trials per condition: ']);
    disp(n_trials_per_cond);
    all_subs_trials_per_cond(i,:)=n_trials_per_cond; %just record for reporting at end
    
    %===========================================%
    % if using ERP data, prepare for classifier %
    %===========================================%
    if erps.use
        erps.tp_inds=mrj_get_tp_inds(erp_trialdata(1).time{1},erps.timepoints_to_use); %assume all conds/trials have same timepoints
        for j=1:n_conds
            %make trial first dim, time second dim, remaining dims can be whatever depending on data type
            data_tmp=zeros([n_trials_per_cond(j) fliplr(size(erp_trialdata(j).trial{1}))]);
            for k=1:n_trials_per_cond(j)
                data_tmp(k,:,:)=(erp_trialdata(j).trial{k})'; %data_tmp indexed by trial, timepoint, channel
            end
            
            if erps.convert_z %convert to z-scores
                for k=1:size(data_tmp,1) %loop thru trials
                    for m=1:size(data_tmp,3) %loop thru channels
                        data_tmp(k,:,m)=mrj_zscores(data_tmp(k,:,m));
                    end
                end
            end
            data(j).erp=mrj_downsample(data_tmp(:,erps.tp_inds,erps.chans_to_use),erps.points_per_bin); %#ok<AGROW>
            %data.erp still indexed by trial, timepoint, channel -- will reshape when it comes time to classify
        end
    end
    
    %======================================================%
    % if using time-frequency data, prepare for classifier %
    %======================================================%
    if freqs.use
        this_ft_cfg=freqs.ft_cfg;  %use this_ft_cfg where old script had mrj_ft_cfg
        %use erp_trialdata where old script had ft_data
        %use n_conds where old script had nconds
        %adapt method to use freqs.chans_to_use instead of nchans
        
        for j=1:n_conds
            this_ft_cfg.toi = min(erp_trialdata(j).time{1}):this_ft_cfg.time_window_delta:max(erp_trialdata(j).time{1});
%             data_tmp=zeros(length(erp_trialdata(j).trial),length(freqs.chans_to_use),length(this_ft_cfg.foi),length(this_ft_cfg.toi));
%             %data_tmp will hold the power spectrum, starts off indexed by trial, channel, frequency, timepoint
            
            this_ft_cfg.channel = erp_trialdata(j).label(freqs.chans_to_use); %make sure this comes out right?
            disp(['Time-frequency calcs, subj ' int2str(i) ', condition ' int2str(j)]); %take out later?
            fa_output = freqanalysis(this_ft_cfg, erp_trialdata(j)); %fa_output.powspec comes out ntrials x nchan x nfreq x ntimepoints
            powspec_alltrials = flipdim(fa_output.powspctrm,3); %still indexed trial,channel,freq,timepoint, but with low frequencies at "bottom"
            is_negtime=fa_output.time<=0;
            powspec_isfinite=isfinite(squeeze(powspec_alltrials(1,1,:,:))); %indexed by freq, timepoint
            if freqs.normalize_to_baseline
                powspec_baseline_inds=false(size(powspec_isfinite));
                for k=1:size(powspec_baseline_inds,1) %loop thru frequencies to get baseline inds for each row
                    baseline_inds=(powspec_isfinite(k,:) & is_negtime);
                    if ~any(baseline_inds)
                        baseline_inds=find(powspec_isfinite(k,:));
                        baseline_inds=baseline_inds(1);
                    end
                    powspec_baseline_inds(k,baseline_inds)=true;
                end
            end
            
            if freqs.normalize_to_baseline && freqs.convert_z
                for k=1:size(powspec_alltrials,1) %loop thru trials
                    for m=1:size(powspec_alltrials,2) %loop thru channels
                        for n=1:size(powspec_alltrials,3) %loop thru frequencies
                            thisrow=squeeze(powspec_alltrials(k,m,n,:));
                            powspec_alltrials(k,m,n,:)=(thisrow - mean(thisrow(powspec_isfinite(n,:)))) / std(thisrow(powspec_isfinite(n,:)));
                            powspec_alltrials(k,m,n,:)=powspec_alltrials(k,m,n,:) - mean(powspec_alltrials(k,m,n,powspec_baseline_inds(n,:)));
                        end
                    end
                end
            elseif freqs.normalize_to_baseline
                for k=1:size(powspec_alltrials,1) %loop thru trials
                    for m=1:size(powspec_alltrials,2) %loop thru channels
                        for n=1:size(powspec_alltrials,3) %loop thru frequencies
                            powspec_alltrials(k,m,n,:)=powspec_alltrials(k,m,n,:) / mean(powspec_alltrials(k,m,n,powspec_baseline_inds(n,:)));
                        end
                    end
                end
            elseif freqs.convert_z
                for k=1:size(powspec_alltrials,1) %loop thru trials
                    for m=1:size(powspec_alltrials,2) %loop thru channels
                        for n=1:size(powspec_alltrials,3) %loop thru frequencies
                            thisrow=squeeze(powspec_alltrials(k,m,n,:));
                            powspec_alltrials(k,m,n,:)=(thisrow - mean(thisrow(powspec_isfinite(n,:)))) / std(thisrow(powspec_isfinite(n,:)));
                        end
                    end
                end
            end
            
            freqs.tp_inds=mrj_get_tp_inds(this_ft_cfg.toi,freqs.timepoints_to_use); %will be different from erps.tp_inds b/c this is post-"resampling" by the time-frequency analysis
            data(j).freq=powspec_alltrials(:,:,:,freqs.tp_inds); %#ok<AGROW>
            %data.freq still indexed by trial,channel,freq,timepoint -- will reshape when it comes time to classify
        end
    end
    
    %==============================================%
    % figure out some preliminary classifier stuff %
    %==============================================%
    n_trials_min=min(n_trials_per_cond);
    n_testtrials_percond=floor(n_trials_min/classifier.kfold);
    n_traintrials_percond=n_testtrials_percond*(classifier.kfold-1);
    traintargs=zeros(n_conds,n_traintrials_percond*n_conds);
    testtargs=zeros(n_conds,n_testtrials_percond*n_conds);
    for j=1:n_conds
        traintargs(j,(j-1)*n_traintrials_percond+1:j*n_traintrials_percond)=1;
        testtargs(j,(j-1)*n_testtrials_percond+1:j*n_testtrials_percond)=1;
    end
    
    %================================%
    % set up data for classification %
    %================================%
    for j=1:n_conds
        if erps.use && ~freqs.use %ERPs only
            data(j).classify=(reshape(data(j).erp,[size(data(j).erp,1) numel(data(j).erp)/size(data(j).erp,1)]))'; %#ok<AGROW> %now nfeatures by ntrials
        elseif ~erps.use && freqs.use %freqs only
            data(j).classify=(reshape(data(j).freq,[size(data(j).freq,1) numel(data(j).freq)/size(data(j).freq,1)]))'; %#ok<AGROW> %now nfeatures by ntrials
        elseif erps.use && freqs.use %both ERPs and freqs
            data(j).classify=[  (reshape(data(j).erp,[size(data(j).erp,1) numel(data(j).erp)/size(data(j).erp,1)]))';
                                (reshape(data(j).freq,[size(data(j).freq,1) numel(data(j).freq)/size(data(j).freq,1)]))';
                             ]; %#ok<AGROW>
        else
            error('Neither ERPs nor frequencies specified for classification!')
        end
        nonfinite_features=any(~isfinite(data(j).classify),2);
        if any(nonfinite_features)
            data(j).classify=data(j).classify(~nonfinite_features,:); %#ok<AGROW>
        end
    end
    
    %===============================%
    % do feature selection, perhaps %
    %===============================%
    if feature_selection.use
        data=mrj_do_feature_selection( data, feature_selection );
    end
    n_features_per_sub(i)=size(data(1).classify,1);

    %==========================%
    % do actual classification %
    %==========================%
    clear allacts allaccs_wta all_wta2_scores; %make sure we start each subject with a clean slate
    for j=1:classifier.nits
        disp(['Subject ' int2str(i) ' of ' int2str(n_subs) ', iteration ' int2str(j) ' of ' int2str(classifier.nits)]);   
        trainpats=zeros(size(data(1).classify,1),n_traintrials_percond*n_conds);
        testpats=zeros(size(data(1).classify,1),n_testtrials_percond*n_conds);
        for k=1:n_conds
            if classifier.shuffle_data_randomly==1 %shuffle data indices within condition
                data(k).classify=data(k).classify(reshape(randperm(numel(data(k).classify)),size(data(k).classify))); %#ok<AGROW>
            elseif classifier.shuffle_data_randomly==2 %just stick in random data the same size as the real data
                data(k).classify=rand(size(data(k).classify)); %#ok<AGROW>
            end
            inds=randperm(size(data(k).classify,2));
            traininds=inds(1:n_traintrials_percond);
            testinds=inds((n_traintrials_percond+1):(n_traintrials_percond+n_testtrials_percond));
            trainpats(:,((k-1)*n_traintrials_percond+1):(k*n_traintrials_percond))=data(k).classify(:,traininds);
            testpats(:,((k-1)*n_testtrials_percond+1):(k*n_testtrials_percond))=data(k).classify(:,testinds);
        end
    
        s=feval(classifier.trainfunc,trainpats,traintargs,classifier.args);
        [acts s]=feval(classifier.testfunc,testpats,testtargs,s);
        if isequal(classifier.testfunc,@test_matlabsvm_mrj)
            allsubs_svm_orig_accs{i,j}=[s(:).accs]; %#ok<AGROW,NASGU>
        end
        allacts(j,:,:)=acts; %#ok<AGROW>
        all_wta2_scores(j,:,:)=tiedrank(acts); %#ok<AGROW>
        acts_maxes=max(acts);
        acts_tmp=zeros(size(acts));
        for k=1:length(acts_maxes)
            acts_tmp(:,k)=(acts(:,k)==acts_maxes(k));
        end
        acts_tmp2=(acts_tmp & testtargs);
        acts_accs_wta=sum(acts_tmp2);
        acts_accs_wta(acts_accs_wta>0)=1./sum(acts_tmp(:,acts_accs_wta>0)); %adjusts for ties -- trials with ties get fractional accuracies
        allaccs_wta(j,:)=acts_accs_wta; %#ok<AGROW>
    end
    
    allsubs_allacts{i}=allacts; %#ok<AGROW,NASGU>
    allsubs_allaccs_wta{i}=allaccs_wta; %#ok<AGROW,NASGU>
    allsubs_all_wta2_scores{i}=all_wta2_scores; %#ok<AGROW,NASGU>
    all_testtargs{i}=testtargs; %#ok<AGROW,NASGU>
end

cd(owd);

keyboard;


%---------------------------------
function tn=mrj_tempname( nchars )

if nargin<1 || isempty(nchars)
    nchars=10;
end

allchars=['a':'z' '0':'9'];

inds=ceil(rand(1,nchars)*length(allchars));
tn=['tmp_' allchars(inds)];

%---------------------------------------------
function tp_inds=mrj_get_tp_inds(tps,tp_range)
    
tp_inds=find(tps>=tp_range(1) & tps<=tp_range(2));

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

%-----------------------------------
function outdata=mrj_zscores(indata)

outdata=(indata - mean(indata))/std(indata);

%---------------------------------------------------------------------
function outdata=mrj_do_feature_selection( indata, feature_selection )

disp(['Doing feature selection: ' feature_selection.type;])
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
