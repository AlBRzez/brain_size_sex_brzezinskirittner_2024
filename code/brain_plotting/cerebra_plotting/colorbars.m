clear
clc

%%
colors_p=readtable("PuOr.csv");
colors_p=table2array(colors_p);
ncolors = colors_p/255;
ncolors2= [[1,1,1];ncolors];
ncolors2(50,:)=[];

%%
%contrasts = ["intercept", "sex_male", "agem", "sex_male_agem"]; %taking vertexwise order for practicity
%range_c = [.7, 1.5, .4, .2]; %volume
%range_c = [.7, 1.5, .3, .1]; %area
%range_c = [.3, .6, .4, .2]; %ct

%contrasts=["intercept", "agem", "sex_male", "sex_male_agem"];%age in months
%range_c = [1, 0.5, 1.5, 0.1]; %used for regular lm trajectories

%contrasts=["intercept", "male_int", "agem", "sex_male", "sex_male_agem", "male_slope"];% allometry - age in months
%range_c=[1.5, 1.5, .1, .15 .15, .1]; %all allometry

%contrasts=["intercept", "agem", "sex_male", "voi", "sex_male_voi"];%age in months; %for fluid inteligence
%range_c = [.1, .1, .25, .25, .1]; %fluid intelligence


% sizes: 478x60; 70x345
%%

for r = [.1, .15, .2, .25, .3, .4, .5, .6, .7, 1, 1.5]
    %r=.25
    figure;imagesc(-r:.01:r);colormap(ncolors);pause(1);set(gcf,'color','w');cb = colorbar('southoutside','Ticks',[-r,0,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
    %figure;imagesc(-r:.01:r);colormap(ncolors);pause(1);set(gcf,'color','w');cb = colorbar('eastoutside','Ticks',[-r,0,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
    colormap(ncolors);pause(1);set(gcf,'color','w');cb = colorbar('eastoutside','Ticks',[-r,0,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
    if r == 1.5
        figure;imagesc(.5:.01:r);colormap(ncolors2);pause(1);set(gcf,'color','w');cb = colorbar('southoutside','Ticks',[.5,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
        %figure;imagesc(.5:.01:r);colormap(ncolors2);pause(1);set(gcf,'color','w');cb = colorbar('eastoutside','Ticks',[.5,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
        colormap(ncolors2);pause(1);set(gcf,'color','w');cb = colorbar('eastoutside','Ticks',[.5,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
    end
end
%%
%close all
%%
r = 1
figure;imagesc(-r:.01:r);colormap(ncolors);pause(1);set(gcf,'color','w');cb = colorbar('southoutside','Ticks',[-r,0,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
colormap(ncolors);pause(1);set(gcf,'color','w');cb = colorbar('eastoutside','Ticks',[-r,0,r],'FontWeight','bold','FontSize',15,'TickLength',0,'Box','off')
