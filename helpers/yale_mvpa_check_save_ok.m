function ok_to_save = yale_mvpa_check_save_ok( fname )

ok_to_save = true;
if exist( fname, 'file' )
    button = questdlg(['The file ' fname ' already exists. Do you want to overwrite?'],'', 'Overwrite', 'No', 'No');
    if strcmp( button, 'Overwrite' )
        ok_to_save = true;
    elseif strcmp( button, 'No' )
        ok_to_save = false;
    else
        error( 'Something wrong in yale_mvpa_check_save_ok.' );
    end
end
