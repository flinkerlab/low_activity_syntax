%% Plot all electrodes from all patients on brain -- project to closest pial surface within anatomical region
clc; clear; close all;

main_path = '/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/';

%% Colors
colors_lightest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/task_gradient_lightest.csv');
colors_lighter = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/task_gradient_lighter.csv');
colors_bright = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/task_line_plots.csv');
colors_darker = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/task_gradient_darker.csv');
colors_darkest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/task_gradient_darkest.csv');

sentence_colors = [
    table2array(colors_bright(strcmp("sp", colors_bright.Var1),2:4));
    table2array(colors_lighter(strcmp("sp", colors_lighter.Var1),2:4));
    table2array(colors_lightest(strcmp("sp", colors_lightest.Var1),2:4));
    1,1,1]
    
list_colors = [1,1,1;
    table2array(colors_lightest(strcmp("lp", colors_lightest.Var1),2:4));
    table2array(colors_lighter(strcmp("lp", colors_lighter.Var1),2:4));
        table2array(colors_bright(strcmp("lp", colors_bright.Var1),2:4))];
cap_d_val = 1;

sentence_cmap = [interp1(linspace(0,cap_d_val, size(list_colors, 1)), list_colors, linspace(cap_d_val,0,128),'pchip');
    interp1(linspace(0,cap_d_val, size(sentence_colors, 1)), sentence_colors, linspace(cap_d_val,0,128),'pchip')];


%% Save colorbar 
% Create a figure
figure;

% Create invisible axes
ax = axes('Visible', 'off');

% Plot dummy data and create colorbar
fake_data = rand(10);
fake_data(1,1) = cap_d_val * 2;
fake_data(1,2) = 0;
fake_data = fake_data - cap_d_val;
imagesc(ax, fake_data); % Plot dummy data
colormap(ax, sentence_cmap); % Set the colormap you want to use
c = colorbar;
c.Visible = 'on'; % Make only the colorbar visible

% Optionally set colorbar properties
set(c, 'Box', 'off');
set(c, 'TickLength', 0);

% Save colorbar
save_cb_dir = fullfile(main_path,"figures","color bar");
if ~exist(save_cb_dir, 'dir')
    mkdir(save_cb_dir);
end

% Save
print(gcf, fullfile(save_cb_dir, strcat("color bar.pdf")), '-dpdf', '-bestfit', '-r1000');

% Close dummy figure
close;


%% Loops
warp_loops = ["warped data", "unwarped data"];
band_loops = ["high_gamma", "beta"];
lock_loops = ["locked_to_production_onset", "locked_to_stimulus_onset"];
bin_loops = ["200ms bins"];%["250ms bins", "200ms bins", "100ms bins"]; % "200ms bins", "150ms bins", "50ms bins"
alpha_loops = ["alpha=0.05"];%["alpha=0.01", "alpha=0.05"]; 
comparison_loops = "sentence_vs_listing"; %, "sentence_vs_naming", "listing_vs_naming"];

for warp_loop_n = 1:length(warp_loops)
    % warp_loop_n = 1;
    warp_loop = warp_loops{warp_loop_n};

    for band_loop_n = 2%:length(band_loops)
        % band_loop_n = 1;
        band_loop = band_loops{band_loop_n};

        for lock_loop_n = 1:length(lock_loops)
            % lock_loop_n = 1;
            lock_loop = lock_loops{lock_loop_n};

            for bin_loop_n = 1:length(bin_loops)
                % bin_loop_n = 1;
                bin_loop = bin_loops{bin_loop_n};

                for alpha_loop_n = 1:length(alpha_loops)
                    % alpha_loop_n = 1;
                    alpha_loop = alpha_loops{alpha_loop_n};
                    
                    for comparison_loop_n = 1:length(comparison_loops)
                        % comparison_loop_n = 1;
                        comparison_loop = comparison_loops{comparison_loop_n};
    
                        disp(strcat("Beginning loops: ",warp_loop," - ",band_loop," - ",lock_loop," - ",bin_loop," - ",alpha_loop,"."))
        
                        % Read in Cohen's Ds and significances
                        current_read_path = fullfile(main_path, "data/", warp_loop, band_loop, "patients - combined", lock_loop, bin_loop);
                        Ds = readtable(fullfile(current_read_path, strcat(comparison_loop,"_cohens_d.csv")));
                        sig = readtable(fullfile(current_read_path, strcat("wilcox_test_",comparison_loop,"_sig_",alpha_loop,".csv")));
                        
                        % Separate elec ID and locations
                        elec_cols = {'elec','region_clinical','MNI_x','MNI_y','MNI_z'};
                        D_elec_info = Ds(:,elec_cols);
                        sig_elec_info = sig(:,elec_cols);
        
                        % Confirm that locations line up between two tables
                        if ~isequal(D_elec_info, sig_elec_info)
                           warning('Electrode locations differ between Cohens D table and significance table (likely different row order).')
                        end
        
                        % Get locations and regions
                        elec_locs = D_elec_info{:,{'MNI_x','MNI_y','MNI_z'}};
                        elec_regions = D_elec_info{:,'region_clinical'};
        
                        % Remove elec cols from data
                        Ds(:,elec_cols) = [];
                        sig(:,elec_cols) = [];
                        
                        % Get windows to loop thru (i.e., remaining column names)
                        window_loops = Ds.Properties.VariableNames;
        
                        % Loop!
                        for window_loop_n = 1:length(window_loops)
                            % window_loop_n = 3;
                            window_loop = window_loops{window_loop_n};
                            
                            % Which elecs sig?
                            current_sig = sig{:,window_loop};
                            current_sig_indices = find(current_sig==1);
    
                            % Subset to just sig
                            current_Ds = Ds{current_sig_indices,window_loop};
                            current_locs = elec_locs(current_sig_indices,:);
                            current_regions = elec_regions(current_sig_indices);
                            
                            %% Get alpha values
                            % Initialize alphaValues with the same size as values
                            alphaValues = zeros(size(current_Ds));
                        
                            % Set alpha to 0 for values in the range [0, 0.2]
                            max_transparent = 0.2; % .2 = negligible D; Cohen 1992
                            alphaValues(abs(current_Ds) <= max_transparent) = 0;
                        
                            % Set alpha to 1 for values 0.5 or higher
                            min_opaque = 0.5; % .5 = medium D; Cohen 1992
                            alphaValues(abs(current_Ds) >= min_opaque) = 1;
                        
                            % Interpolate alpha linearly for values between 0.2 and 0.5
                            % Calculate the slope for linear interpolation (rise/run)
                            slope = 1 / (min_opaque - max_transparent);
                            % Apply linear interpolation
                            for i = 1:length(current_Ds)
                                if abs(current_Ds(i)) > max_transparent && abs(current_Ds(i)) < min_opaque
                                    alphaValues(i) = slope * (abs(current_Ds(i)) - max_transparent);
                                end
                            end
                            
    
                            %% Get color values
                            colorValues = NaN(size(current_Ds, 1), 3);
                            for i = 1:length(current_Ds)
                                % i = 1
                                [~, current_color] = min(abs(linspace(-cap_d_val, cap_d_val, size(sentence_cmap, 1))' - current_Ds(i)));
                                colorValues(i,:) = sentence_cmap(current_color,:);
                            end
                            
                            %% Plot!
                            cd '/Users/am4611/Dropbox/Research/ChickenSyntax/brain plot code from amir/visualization-tools-main/matlab/'
                            %AnnotFile='./SampleData/MNI/ch2_template.lh.aparc.split_STG_MTG.annot';
                            %BrainFile='./SampleData/MNI/ch2_template_lh_pial_120519.mat';
                            AnnotFile='./SampleData/MNI-FS/FSL_MNI152.lh.aparc.split_STG_MTG.annot';
                            BrainFile='./SampleData/MNI-FS/FSL_MNI152_lh_pial.mat';
                            
                            % Project electrodes to nearest pial surface within region
                            [elec_locs_proj, ~] = project_elecs_pial(current_locs, current_regions, BrainFile, AnnotFile);
                            
                            % Initialize brain
                            VT = visualtools('Subj', 'MNI',...
                                             'HS', 'lh',...
                                             'flag_UseAnnots', false,...
                                             'BrainFile', BrainFile,...
                                             'AnnotFile', AnnotFile);
        
                            %% Plot all elecs, colored by whether active
                            % plot electrodes on the brain
                            p=VT.PlotElecOnBrain(elec_locs_proj, ...
                                               'ElecColor', colorValues,...
                                               'ElecAlpha', ones(size(alphaValues)),...
                                               'BrainColor', [1,1,1],...
                                               'cmap',sentence_cmap,...
                                               'clim',[0,cap_d_val],...
                                               'radius',2.36);
                            p.AmbientStrength =0.5; %def 0.3
                            p.DiffuseStrength =0.6; % def 0.8
                            %p.Visible = 'off'
                            colorbar('off')
    
                            % Set aspect ratio of PDF
                            ax = gca; % Get current axes
                            aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
                            figureWidth = 8; % Define the desired width of the figure in inches
                            figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
                            set(gcf, 'PaperUnits', 'inches'); % Set figure size
                            set(gcf, 'PaperSize', [figureWidth figureHeight]);
                            
                            % Set up save directory
                            save_fig_dir = fullfile(main_path,"figures",warp_loop,band_loop,lock_loop,bin_loop,alpha_loop);
                            if ~exist(save_fig_dir, 'dir')
                                mkdir(save_fig_dir);
                            end
                            
                            % Save
                            print(gcf, fullfile(save_fig_dir, strcat("window_loop - ",window_loop," - lower cap.pdf")), '-dpdf', '-bestfit', '-r1000');
                            
                            close;
    
                        end
                    end
                end
            end
        end
    end
end

