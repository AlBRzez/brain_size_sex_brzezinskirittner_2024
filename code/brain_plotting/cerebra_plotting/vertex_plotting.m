clear
clc
%% Loading tools
tools = '/data/zeiyas/tools/';
addpath(genpath(strcat(tools, 'matlab_toolboxes/export_fig-master/')))
addpath(genpath(strcat(tools,'matlab_toolboxes/surfstat')))
addpath(genpath(strcat(tools,'useful_scripts')))

% read surface template
model_folder_s = strcat(tools,'parcellations_atlas_mni/models/');
s = SurfStatReadSurf( {strcat(model_folder_s,'surface_models/icbm_avg_mid_sym_mc_left.obj'),...
    strcat(model_folder_s,'surface_models/icbm_avg_mid_sym_mc_right.obj')} );
% make Cerebra surface and fill in the holes
atlases = load(strcat(tools,'parcellations_atlas_mni/yz_atlasses/all_atlases_vectorized.mat'));

%%
colors_p=readtable("PuOr.csv");
colors_p=table2array(colors_p);

%% Loading data
folder='/data/zeiyas/brzali/brain_size_sex_brzezinskirittner_2024/outputs/plots/brainplots';
Info=readtable(strcat(folder, '/../../vertexwise/b_volume_civet_all_20'),'Delimiter',',');%for


%%
out_fold='vertex_vol';
model="vol";

% Vertexwise color ranges
range_c = [.7, 1.5, .4, .2]; %volume
%range_c = [.7, 1.5, .3, .1]; %area
%range_c = [.3, .6, .4, .2]; %ct

%% regions matrix
% Each column represents the regional estimates for one sample for one contrast
samples=["extreme", "random", "age_mat", "matched"];
n_samp=length(samples);
contrasts = ["intercept", "sex_male", "agem", "sex_male_agem"]; %taking vertexwise order for practicity
n_cont=length(contrasts);

%% plotting in a loop
ncolors = colors_p/255;
for i=1:n_samp
    for j=1:n_cont
        out = Info(:,(i-1)*n_cont+j);
        out_s = table2array(out);
        out_s(abs(out_s) < range_c(j)/93 & out_s ~= 0) = range_c(j)/93 + eps;

        %if i == 4
        %    colormap(ncolors)
        %    figure;SurfStatViewData_yz_22(out_s,s,[-range_c(j),range_c(j)],strcat('estimate'));
        %    pause(1);set(gcf,'color','w')
        %    colormap(ncolors)
        %else
        colormap(ncolors)
        figure;SurfStatViewData_yz_22(out_s,s,[-range_c(j),range_c(j)],strcat('estimate'));
        pause(1);set(gcf,'color','w');colorbar off
        colormap(ncolors)
        %end

        export_fig(strcat(folder, '/', out_fold, '/', model, '_', samples(i),'_',contrasts(j),'.png'),'-m4')
    end
end

%%
%figure;imagesc(-1:.01:1);colormap(nrdbu);pause(1);set(gcf,'color','w');colorbar
%export_fig(strcat(folder, '/colorbar.png'))
close all