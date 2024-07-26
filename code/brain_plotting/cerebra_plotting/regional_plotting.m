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

%% Loading data
folder='/data/zeiyas/brzali/brain_size_sex_brzezinskirittner_2024/outputs/plots/brainplots';
Info=readtable(strcat(folder, '/brainplots_data/allom_60m.csv'),'Delimiter',',');%for
CerebraInfo=readtable('cerebra_reference.csv','Delimiter',',');CerebraInfo = sortrows(CerebraInfo,'newlabels','ascend');
%%
%rdbu=readtable("redblue_rgb.csv");
rdbu=readtable("redblue_rgb.csv");
rdbu=table2array(rdbu);

%%
out_fold='allom_60m_reg_colorbar';

%% regions matrix
% Each column represents the regional estimates for one sample for one contrast
samples=["extreme", "random", "age_mat", "matched"];
%samples=["random", "matched"]; % for fluid inteligence
n_samp=length(samples);
%contrasts=["intercept", "agem", "sex_male", "sex_male_agem"];%age in months
contrasts=["intercept", "male_int", "agem", "sex_male", "sex_male_agem", "male_slope"];% allometry - age in months
contrasts=unique(Info.term_c); %for fluid inteligence
n_cont=length(contrasts);
%unique_models=unique(Info.mod);
model="lin_int_regular";

%% For the cool iterations
% looping through each sample and contrast to get the ordered values
n=1;
for smp=1:n_samp
    samp=Info(ismember(Info.df, [samples(smp)]) & ismember(Info.mod, model), :);
    for i=1:n_cont
        tmp=samp(ismember(samp.term_c, [contrasts(i)]), :);
        tmp = sortrows(tmp,'newlabels','ascend');
        regions(:,n)=table2array(tmp(:,"estimate"));
        n=n+1;
    end
end

%% surface visualization
Cerebra_s = ICBM_volume_to_surface_map(atlases.all_atlases_volume.Cerebra_gm_h);

Cerebra_s_stat = tabulate(Cerebra_s);
remove_threshold = 10;remove_regions = Cerebra_s_stat((Cerebra_s_stat(:,2)<remove_threshold),1);
Cerebra_s(ismember(Cerebra_s,remove_regions))=0;

complementary_atlas= atlases.all_atlases_surface.Schaefer_7_1000;
Cerebra_s(complementary_atlas==0)=0;
ind_miss = find(Cerebra_s==0 & complementary_atlas~=0);

for i=1:size(ind_miss)
    temp_cerebra = Cerebra_s(complementary_atlas==complementary_atlas(ind_miss(i)));
    temp_val(i,1)= mode(temp_cerebra(temp_cerebra~=0));
end
Cerebra_s(ind_miss)=temp_val;clear temp_val
Cerebra_s(isnan(Cerebra_s))=0;
%% project results to the surface cerebra
u_regions = unique(Cerebra_s);u_regions(u_regions==0)=[];

%% plotting in a loop
% Regional color ranges for each analysis
%range_c = [1, 0.5, 1.5, 0.1]; %used for regular lm trajectories
range_c=[1.5, 1.5, .1, .15 .15, .1]; %all allometry
%range_c = [.1, .1, .2, .1, .25]; %fluid intelligence

nrdbu = rdbu/255;
nrdbu2 = [[1,1,1];nrdbu];
nrdbu2(48,:)=[];
%%
for i=1:n_samp
    for j=1:n_cont
        out = regions(:,(i-1)*n_cont+j);out = out(u_regions);
        out_s = region_to_atlas(out,Cerebra_s);
        out_s(abs(out_s) < range_c(j)/93 & out_s ~= 0) = range_c(j)/93 + eps;
        %cmap=colormap('jet');cmap_new = [cmap(1:32,:);[192,192,192]/255;cmap(33:64,:)]; colormap(cmap_new)
        %if ismember(j, [1, 2])
        %    colormap(nrdbu2)
        %    figure;SurfStatViewData_yz_22(out_s,s,[.5,range_c(j)],strcat('estimate'));
        %    pause(1);set(gcf,'color','w');colorbar off
        %    colormap(nrdbu2)
        %else
        if ismember(j, [1,2])
            cm = nrdbu2;
        else
            cm = nrdbu;
        end

        if i == 4
            colormap(cm)
            figure;SurfStatViewData_yz_22(out_s,s,[-range_c(j),range_c(j)],strcat('estimate'));
            pause(1);set(gcf,'color','w')
            colormap(cm)
        else
            colormap(cm)
            figure;SurfStatViewData_yz_22(out_s,s,[-range_c(j),range_c(j)],strcat('estimate'));
            pause(1);set(gcf,'color','w');colorbar off
            colormap(cm)
        end
        %end


        export_fig(strcat(folder, '/', out_fold, '/', model, '_', samples(i),'_',contrasts(j),'.png'),'-m4')
    end
end
%%
%imagesc(-1:.01:1);colormap(nrdbu);colorbar