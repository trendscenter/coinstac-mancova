%function [] = remote_mat_script(stats_type_arr, input_dir, out_dir)
%% TODO: Determine what to return
% ------------------------------------------------
% At master site, load new mat files and perform
% ------------------------------------------------
fprintf('Inside remote matlab script \n')
%%% Docker test
%input_dir='/computation/temp_output/';
%output_dir='/computation/temp_output/aggregated_mancovan/'


%% Local mac test
%input_dir='/Users/sbasodi1/Downloads/mancova_srinivas/spectra_test/';
%output_dir='/Users/sbasodi1/Downloads/mancova_srinivas/spectra_test/merged';

input_dir='/input/remote/simulatorRun/'
output_dir='/output/remote/simulatorRun/coinstac_reduced_mancovan_univariate/'

if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end


%stats_type='SPECTRA';
stats_type_arr= {'FNC', 'SPECTRA'};
allfiles = dir(fullfile(input_dir,'local*/*_reduced_stats_mancovan*.mat'));
num_local_clients=numel(allfiles);
all_terms='';
regressors='';
comp_network_names={};

fnc_stats_arr=[];
spectra_stats_arr=[];


fprintf('Num files loaded: %i \n',numel(allfiles))

%Assert len of all the arrays are the same
for fileid = 1:numel(allfiles)
    fprintf('Filename: %s\n',fullfile(input_dir, allfiles(fileid).name))
    if fileid == 1
        temp_mat=load(fullfile(allfiles(fileid).folder, allfiles(fileid).name));
        spectra_dims=size(temp_mat.spectraStatsArr);
        spectra_stats_arr=cell(num_local_clients, spectra_dims(1), spectra_dims(2));
        all_terms=temp_mat.allTerms;
        regressors=temp_mat.regressors;
        for col = 1:size(temp_mat.compNetworkNames,2)
            comp_network_names(col,:)={temp_mat.compNetworkNames(col).name, temp_mat.compNetworkNames(col).value};
        end
    end
    if any(strcmpi(stats_type_arr,'FNC'))
        fnc_stats_arr = [fnc_stats_arr; temp_mat.fncStatsArr];
    end
        
    if any(strcmpi(stats_type_arr,'SPECTRA'))
        spectra_stats_arr(fileid, :, :) = temp_mat.spectraStatsArr;
    end
        
    if any(strcmpi(stats_type_arr,'SM'))
        %TODO: Deal with SM stats; Has 53 mat files corresponding to each fnc comp
        disp('TODO: Not implemented SM stats');
    end
end


%if any(strcmp(stats_type_arr,'FNC'))
%    fprintf('\nProcessing STATS TYPE: FNC\n')
%    plot_FNC(fnc_stats_arr, all_terms, comp_network_names, regressors, output_dir)
%end

%if any(strcmp(stats_type_arr,'SPECTRA'))
%    fprintf('\nProcessing STATS TYPE: SPECTRA\n')
%    plot_Spectra(spectra_stats_arr, all_terms, comp_network_names, regressors, output_dir)
%end

%if any(strcmp(stats_type_arr,'SM'))
%    %TODO: Deal with SM stats; Has 53 mat files corresponding to each fnc comp
%    a='TODO: Not implemented'
%end

%function plot_Spectra(spectra_stats_arr, all_terms, comp_network_names, regressors, output_dir)
if any(strcmp(stats_type_arr,'SPECTRA'))
    load icatb_colors coldhot;
    [num_clients,num_comp, num_terms]=size(spectra_stats_arr)
    
    for curr_term = 1:num_terms
        
        t_r_allcomp=[];
        p_r_allcomp=[];
        stats_r_allcomp=[];

        for curr_comp = 1:num_comp
            curr_stats_arr=spectra_stats_arr(:,curr_comp,curr_term);
            termNo=curr_stats_arr{1,1}.Term

            %TODO check with Srinivas why termNo [1,2] errors
            if length(curr_term) == 1
                    %fprintf('Processing Term: %i\n', termNo)
                    %TODO check with Srinivas why termNo [1,2] errors
                    if length(termNo) == 1
                        [t_r, p_r, stats_r] = mT_reduce(curr_stats_arr, all_terms, termNo, {'verbose'});
                        new_t_r=-sign(t_r).*log10(p_r + eps);
                        t_r_allcomp=[t_r_allcomp;new_t_r];
                        p_r_allcomp=[p_r_allcomp;p_r];
                        stats_r_allcomp=[stats_r_allcomp;stats_r];
                    end
            end
        end

        if ~ isempty(t_r_allcomp)
            fprintf('Finished mT_reduce on all components; Generating images \n');
            F=figure('visible','off'); 
            imagesc(t_r_allcomp);colormap(coldhot);
            title_str=sprintf('Spectra Component Map for Term: %i (%s)',termNo, string(regressors(termNo)));
            title(title_str);
            xlabel('t_r vector');
            ylabel('Components');

            output_file_name=fullfile(output_dir,sprintf('NEW_SPECTRA_modtr_corr_term_%i_%s.png',termNo, string(regressors(termNo))));

            saveas(F,output_file_name);
            clf(F);
            close(F);
        end
    end
end

%function plot_FNC(fnc_stats_arr, all_terms, comp_network_names, regressors, output_dir)
if any(strcmp(stats_type_arr,'FNC'))

    load icatb_colors coldhot;
    [num_clients,num_terms]=size(fnc_stats_arr)

    for col = 1:num_terms
        curr_stats_arr=fnc_stats_arr(:,col);
        termNo=curr_stats_arr{1,1}.Term
        fprintf('Processing Term: %i\n', termNo)
        %TODO check with Srinivas why termNo [1,2] errors
        if length(termNo) == 1
            [t_r, p_r, stats_r] = mT_reduce(curr_stats_arr, all_terms, termNo, {'verbose'});
            fprintf('Finished mT_reduce; Generating images \n')

            %%IMAGE GENERATION CODE FOR FNC
            % network names
            network_names = comp_network_names(:, 1);

            % network Lengths
            comps = [comp_network_names{:, 2}];
            comps = comps(:)';
            network_values = cellfun(@length, comp_network_names(:, 2), 'UniformOutput', false);
            network_values = [network_values{:}];

            new_t_r=-sign(t_r).*log10(p_r + eps);
            M = icatb_vec2mat(new_t_r);

            % FNC matrix (53 x 53)
            %M = icatb_vec2mat(t_r);

            % axesTitle = 'FNC correlations';
            axes_title = sprintf('FNC correlations Term:%i - %s',termNo, string(regressors(termNo)));

            CLIM = max(abs(M(:)));
            CLIM = 1
            CLIM = -log10(eps)

            gH = icatb_getGraphics(axes_title, 'graphics',  'FNC Correlations', 'on');
            set(gH, 'resize', 'on');
            axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
            [F,A,C,I] = icatb_plot_FNC(M, [-CLIM, CLIM], cellstr(num2str(comps(:))), (1:length(comps)), gH, axes_title, axesH, ...
                network_values, network_names);
            %colormap(jet(64));
            colormap(coldhot);
            output_file_name=fullfile(output_dir,sprintf('NEW_FNC_modtr_corr_term_%i_%s.png',termNo, string(regressors(termNo))))

            saveas(F,output_file_name);
            clf(F);
            close(F);
        end

    end
    
end

