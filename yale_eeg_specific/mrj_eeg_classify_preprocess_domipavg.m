function erp_data = mrj_eeg_classify_preprocess_domipavg( paths, erps, subnum )

owd=pwd;

%--------------------------------------
% SET UP DATA / PARAMETERS FOR MIPAVG
%--------------------------------------

% get name of folder with .EEG files in it, and then the files themselves
eeg_dir = fullfile( paths.experiment_directory, paths.sub_list{subnum}, paths.eeg_subdir );
cd(eeg_dir);
eeg_files=cellstr(spm_select('List',pwd,'.EEG$'));

% high-pass filter or detrending?
flags={'-triponly'};                                        % assume we don't care about saving averaged ERPs
if ischar( erps.highpass_filter ) && strcmpi( erps.highpass_filter, 'detrend' )
    flags{end+1}='-lin';
elseif isnumeric ( erps.highpass_filter ) && erps.highpass_filter > 0
    flags{end+1}=['-f' num2str( erps.highpass_filter ) ',0'];
elseif ~isempty( erps.highpass_filter ) && ~( isnumeric( erps.highpass_filter ) && erps.highpass_filter == 0 )
    error('Unrecognized value for erps.highpass_filter');
end

% generate a suitable temporary filename for fieldtrip-style input
mip_output_tempname=mrj_tempname;
while exist([mip_output_tempname '.avg'],'file') || exist([mip_output_tempname '.log'],'file') || exist([mip_output_tempname '.hdr'],'file') || exist([mip_output_tempname '_fieldtrip.mat'],'file')
    mip_output_tempname=mrj_tempname;
end

% run mipavg to pre-process trials and split into epochs
mrj_mipavg3_inds( eeg_files{:}, paths.p_arf, paths.p_sgc, paths.p_mon, [mip_output_tempname '.avg'], [mip_output_tempname '.log'], flags{:} );

% load output data and delete temporary file
erp_data = load([mip_output_tempname '_fieldtrip.mat']);
erp_data = erp_data.mipavg_fieldtrip_data;
delete( [mip_output_tempname '_fieldtrip.mat'] );

cd(owd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTILITY FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------
%---------------------------------
function tn=mrj_tempname( nchars )

if nargin<1 || isempty(nchars)
    nchars=10;
end

allchars=['a':'z' '0':'9'];

inds=ceil(rand(1,nchars)*length(allchars));
tn=['tmp_' allchars(inds)];
