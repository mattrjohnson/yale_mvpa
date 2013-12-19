%results for regular WTA
for zzz=1:length(allsubs_allaccs_wta)
yyy(zzz)=mean(mean(allsubs_allaccs_wta{zzz}));
end
yyy'
clear yyy zzz
%results for WTA2
for zzz=1:length(allsubs_all_wta2_scores)
tmp=squeeze(mean(allsubs_all_wta2_scores{zzz},1));
yyy(zzz)=mean(tmp(logical(all_testtargs{zzz})));
end
yyy'/10
clear tmp yyy zzz;
%results for SVM raw accuracies
for zzz=1:size(allsubs_svm_orig_accs,1)
for yyy=1:size(allsubs_svm_orig_accs,2)
xxx(zzz,yyy)=mean(mean(allsubs_svm_orig_accs{zzz,yyy}));
end
end
mean(xxx,2)/2
clear xxx yyy zzz