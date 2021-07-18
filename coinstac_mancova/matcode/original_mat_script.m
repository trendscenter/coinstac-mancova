addpath(genpath('/Users/sbasodi1/Downloads/mancova_srinivas/gift-master/'))

baseDir='/Users/sbasodi1/workplace/mancova/mancova_suni_test_mac/success_output_dirs/server_run3_allsub_fnc_spectra/brad_orig_output/'
root_regression_dir=fullfile(baseDir,'output/remote/simulatorRun/coinstac-univariate-regression/')
%root_regression_dir=fullfile(baseDir,'output/remote/simulatorRun/coinstac-univariate-ttest2-gender/')
output_dir=fullfile(baseDir,'aggregated_mancovan/')
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end
mancovan_mat= load(fullfile(root_regression_dir,'coinstac-gica_merge_mancovan.mat'));
regressors=mancovan_mat.mancovanInfo.regressors;
allTerms=mancovan_mat.mancovanInfo.terms;
compNetworkNames=mancovan_mat.mancovanInfo.comp;

%process_FNC_stats(root_regression_dir, allTerms, compNetworkNames, regressors, output_dir)
process_SPECTRA_stats(root_regression_dir, allTerms, compNetworkNames, regressors, output_dir)



function process_SPECTRA_stats(root_regression_dir, allTerms, compNetworkNames, regressors, output_dir)
    load icatb_colors coldhot;
    
    spectra_files=dir(fullfile(root_regression_dir, 'spectra_stats/coinstac-gica_merge_mancovan_results_spectra_*.mat'));
    %Assert len of all the arrays are the same from input: numel(spectra_files)
    spectra_stats_arr=[];
    spectra_data_arr=[];
    for fileid = 1:numel(spectra_files)
        temp_mat=load(fullfile(spectra_files(fileid).folder, spectra_files(fileid).name));
        if fileid == 1
            spectra_data_dims=size(temp_mat.spectra_tc);
            spectra_data_arr=cell(numel(spectra_files), spectra_data_dims(1), spectra_data_dims(2));
        end
        
        spectra_stats_arr=[spectra_stats_arr; temp_mat.UNI.stats];
        %spectra_data_arr(fileid, :, :)= temp_mat.spectra_tc;
    end
        
    [num_comp, num_terms]=size(spectra_stats_arr)
    
    for curr_term = 1:num_terms
        t_r_allcomp=[];
        p_r_allcomp=[];
        stats_r_allcomp=[];

        for curr_comp = 1:num_comp
            curr_stats_arr=spectra_stats_arr(curr_comp,curr_term);
            termNo=curr_stats_arr{1,1}.Term

            %TODO check with Srinivas why termNo [1,2] errors
            if length(curr_term) == 1
                    %fprintf('Processing Term: %i\n', termNo)
                    %TODO check with Srinivas why termNo [1,2] errors
                    if length(termNo) == 1
                        %[t_r, p_r, stats_r] = mT_reduce(curr_stats_arr, all_terms, termNo, {'verbose'});
                        temp_mat=load(fullfile(spectra_files(curr_comp).folder, spectra_files(curr_comp).name));
                        data= temp_mat.spectra_tc;
                        
                        X_reduced=curr_stats_arr{1,1}.X;
                        
                        [t_r, p_r, stats_r] = mT(data, X_reduced, allTerms, termNo, {'verbose'});
                        new_t_r=-sign(t_r).*log10(p_r + eps);
                        new_t_r=t_r;
                        t_r_allcomp=[t_r_allcomp;new_t_r];
                        p_r_allcomp=[p_r_allcomp;p_r];
                        stats_r_allcomp=[stats_r_allcomp;stats_r];
                    end
            end
        end

        if ~ isempty(t_r_allcomp)
            fprintf('Finished mT on all components; Generating images \n');
            F=figure('visible','off'); 
            imagesc(t_r_allcomp);colormap(coldhot);
            title_str=sprintf('Spectra Component Map for Term: %i (%s)',termNo, string(regressors(termNo)));
            title(title_str);
            xlabel('t_r vector');
            ylabel('Components');

            output_file_name=fullfile(output_dir,sprintf('OLD_SPECTRA_origtr_corr_term_%i_%s.png',termNo, string(regressors(termNo))))

            saveas(F,output_file_name);
            clf(F);
            close(F);
        end
    end
    
    
end

function process_FNC_stats(root_regression_dir, allTerms, compNetworkNames, regressors, output_dir)
    load icatb_colors coldhot;

    fnc_stats_mat= load(fullfile(root_regression_dir,'fnc_stats/coinstac-gica_merge_mancovan_results_fnc.mat'));
    %fncStatsArr=fnc_stats_mat.UNI.stats;
    data=fnc_stats_mat.fnc_corrs;
    %%%% CODE TO REPEAT FOR EACH TERM
    for statsnum = 1:numel(fnc_stats_mat.UNI.stats)
        X_reduced=fnc_stats_mat.UNI.stats{1,statsnum}.X;
        termNo=fnc_stats_mat.UNI.stats{1,statsnum}.Term

        [ t_u, p_u, stats_u] = mT(data, X_reduced, allTerms, termNo, {'verbose'});
        for col = 1:size(compNetworkNames,2)
          comp_network_names(col,:)={compNetworkNames(col).name, compNetworkNames(col).value};
        end
        network_names = comp_network_names(:, 1);
        % network Lengths
        comps = [comp_network_names{:, 2}];
        comps = comps(:)';
        network_values = cellfun(@length, comp_network_names(:, 2), 'UniformOutput', false);
        network_values = [network_values{:}];

        new_t_u=-sign(t_u).*log10(p_u + eps);
        new_t_u=t_u
        M = icatb_vec2mat(new_t_u);
        CLIM = max(abs(M(:)));
        CLIM = -log10(eps)
        CLIM = 1

        axesTitle = sprintf('FNC correlations Term:%i - %s',termNo, string(regressors(termNo)));
        gH = icatb_getGraphics(axesTitle, 'graphics',  'FNC Correlations', 'on');
        set(gH, 'resize', 'on');
        axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        [F,A,C,I] = icatb_plot_FNC(M, [-CLIM, CLIM], cellstr(num2str(comps(:))), (1:length(comps)), gH, axesTitle, axesH, ...
        network_values, network_names);
        %colormap(jet(64));
        colormap(coldhot);
        
        output_file_name=fullfile(output_dir,sprintf('OLD_FNC_origtr_CLIM_1_corr_term_%i_%s.png',termNo, string(regressors(termNo))))
        saveas(F,output_file_name);
        clf(F);
        close(F);
    end
end

