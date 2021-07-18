%mat_info=load('/computation/mancova_structure_info.mat');
root_regression_dir=mat_info.univariateResultsDir;
file_prefix=mat_info.mancovanInputFilesPrefix

%try
%    mancovan_mat= load(fullfile(root_regression_dir,'coinstac-gica_mancovan.mat'));
%    file_prefix='coinstac-gica_';
%catch ME
%    mancovan_mat= load(fullfile(root_regression_dir,'coinstac-gica_merge_mancovan.mat'));
%    file_prefix='coinstac-gica_merge_';
%end

mancovan_mat= load(fullfile(root_regression_dir,strcat(file_prefix,'mancovan.mat')));
regressors=mancovan_mat.mancovanInfo.regressors;
allTerms=mancovan_mat.mancovanInfo.terms;
compNetworkNames=mancovan_mat.mancovanInfo.comp;

for f_indx = 1 : length(mat_info.featureStats)
    curr_feature= mat_info.featureStats{f_indx};
    if strcmpi(curr_feature, 'fnc correlations')
        %%%%%%%%%%FNC stats
        disp('Adding FNC stats ..');

        fnc_stats_mat= load(fullfile(root_regression_dir, strcat('fnc_stats/',file_prefix,'mancovan_results_fnc.mat')));
        fncStatsArr=fnc_stats_mat.UNI.stats;

    elseif (strcmpi(curr_feature, 'timecourses spectra'))
        disp('Adding spectra stats ..');
        %%%%%%%%%%Spectra Stats
        %spectra_files=dir(fullfile(root_regression_dir, 'spectra_stats/coinstac-gica_mancovan_results_spectra_*.mat'));
        spectra_files=dir(fullfile(root_regression_dir,  strcat('spectra_stats/', file_prefix, 'mancovan_results_spectra_*.mat')));

        spectraStatsArr=[];
        for fileid = 1:numel(spectra_files)
            temp_mat=load(fullfile(spectra_files(fileid).folder, spectra_files(fileid).name));
            stats_arr=temp_mat.UNI.stats;
            spectraStatsArr=[spectraStatsArr; temp_mat.UNI.stats];
        end
    end
end


%TODO: Later add SM stats

output_file_name=strcat(mat_info.clientID,'_reduced_stats_mancovan');
%% Save the extracted stats/info to output dirs
for dir_indx = 1 : length(mat_info.outputDirs)
    output_dir=mat_info.outputDirs{dir_indx};
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    output_file=fullfile(output_dir,output_file_name);
    disp(sprintf('saving output in: %s', output_file));
    save(output_file, 'allTerms', 'regressors', 'compNetworkNames', 'fncStatsArr', 'spectraStatsArr');

end

