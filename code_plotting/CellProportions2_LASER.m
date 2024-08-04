%% Cell properties plots (asymmetry, shifts) - NO VELOCITY CONTROL
close all
clear


direction_name = {'UP','DOWN'};
laser_name = {'OFF','ON'};

cell_name = {'LOCKING','PRECESSING'};

% data_origin = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\!LaserTrondheim\Matteo_GLMLaser_speed\';
 data_origin = 'F:\!papers\!LaserTrack_Paper\!DATA_LaserTrondheim\CorrectData\';
 %control
%  data_origin = 'F:\processing\!Laser_HD\control_animals\CTRL_Data_2\'

% prop_dir = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\';
% glm_dir = 'C:\Users\Matteo\Desktop\LaserTrondheim\Matteo_GLMLaser\';


%Chose which cells to take depending of remap of PEAK 
remap_value = 15;
gro = 0;

    
theta_All_OFF = [];
theta_All_ON = [];
pf_l_all = [];
pfsasy_All = [];

for Laser = [0 1]

Sparsity_All = [];
Skaggs_All = [];
pfs_All = [];
MeanR_All = [];
MaxR_All = [];
MPos_All = [];


for anm = [1 2 3 4 5 6]
for sess = 1:3
            
   if anm == 1 && sess == 2 || anm == 1 && sess == 3 || anm == 2 && sess == 2 || anm == 2 && sess == 3 || anm == 3 && sess == 3 || anm == 4 && sess == 1 || anm == 5 && sess == 1 || anm == 11 && sess == 3 || anm == 12 && sess == 2 || anm == 12 && sess == 3 || anm == 13 && sess == 3
    continue   
   end 
    
for direction = 1:2

    
     disp(['Doing Something_' direction_name{direction} '_' laser_name{Laser+1} '_A' num2str(anm) '_S' num2str(sess) '_F' num2str(1)])
    
load([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) 'Sp2'],...
    'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')    

% load(['/Users/federico/GitHub/Basins_GLM_paper/Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],...
%     'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')   
    
Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 

if(direction == 2)
   pf_asym = -pf_asym; 
   pf_limits = 100 - fliplr(pf_limits);
   
   
end

pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>0 & pf_limits(:,2)<100)); 

if Laser ==0 
    theta_All_OFF = cat(1,theta_All_OFF,theta_scores);
elseif Laser == 1
    theta_All_ON = cat(1,theta_All_ON,theta_scores);
end 

pf_l_all = cat(1,pf_l_all,pf_limits);

size(theta_scores)
size(MaxRate)


end
end
end 

end


n_cells = size(pf_l_all,1)/2;
pf_1 = pf_l_all(1:n_cells,:);
pf_2 = pf_l_all(n_cells+1:end,:);

pa_1 = pfsasy_All(1:n_cells);
pa_2 = pfsasy_All(n_cells+1:end);


for cond = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

if(cond==1)
sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cond==2)
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end

figure(301)
subplot(2,1,cond)
histogram(pa_1(sel),-10:0.5:10,'Normalization','probability')
hold on
histogram(pa_2(sel),-10:0.5:10,'Normalization','probability')
xlabel('Asymmetry')
ylabel('Probability')
title(cell_name{cond})


linkaxes


figure(302)
subplot(1,2,cond)
[F,X]=ecdf(pa_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pa_2(sel));

plot(X,F,'LineWidth',2)
xlabel('Asymmetry')
ylabel('Cumulative Probability')
[h(cond), p(cond)] = ttest(pa_1(sel),pa_2(sel))
title([cell_name{cond} ': h=' num2str(h(cond)) ': p=' num2str(p(cond))])
linkaxes


figure(5002)
subplot(1,2,cond)
[F,X]=ecdf(pf_2(sel,1)-pf_1(sel,1));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pf_2(sel,2)-pf_1(sel,2));
plot(X,F,'LineWidth',2)

[F,X]=ecdf(pf_2(sel,3)-pf_1(sel,3));
plot(X,F,'LineWidth',2)
xlim([-50 50])
plot([0 0],[0 1],'--k')
plot([-50 50],[0.5 0.5],'--k')

xlabel('Shift (cm)')
ylabel('Cumulative Probability')
title(cell_name{cond})
legend({'Start','Centre','End'},'FontSize',16,'Location','northwest')

end

figure(301)
legend({laser_name{1}, laser_name{2}},'FontSize',16)

figure(302)
legend({laser_name{1}, laser_name{2}},'FontSize',16,'Location','east')


for cond = 1


if(cond==1)
    cltype = find(theta_All_OFF<0 );
elseif(cond==2)
    cltype = find(theta_All_OFF>0);
end

figure(5001)
subplot(1,3,1)
histogram(pf_2(cltype,1)-pf_1(cltype,1),-50:4:50,'Normalization','probability')
mean_start = mean((pf_2(cltype,1)-pf_1(cltype,1)))
pf_start = (pf_2(cltype,1)-pf_1(cltype,1));
hold on
ms = plot(mean_start,0);
set(ms,'Marker','square',...
    'MarkerFaceColor',[1 .6 .6])
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field Start')

subplot(1,3,2)
histogram(pf_2(cltype,2)-pf_1(cltype,2),-50:4:50,'Normalization','probability')
mean_center = mean((pf_2(cltype,2)-pf_1(cltype,2)))
pf_center = (pf_2(cltype,2)-pf_1(cltype,2));
hold on
mc = plot(mean_center,0);
set(mc,'Marker','square',...
    'MarkerFaceColor',[1 .6 .6])
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field Centre')

subplot(1,3,3)
histogram(pf_2(cltype,3)-pf_1(cltype,3),-50:4:50,'Normalization','probability')
mean_end = mean((pf_2(cltype,3)-pf_1(cltype,3)))
pf_end = (pf_2(cltype,3)-pf_1(cltype,3));
hold on
me = plot(mean_end,0);
set(me,'Marker','square',...
    'MarkerFaceColor',[1 .6 .6])
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field End')

linkaxes
end 

figure(6001)
subplot(1,2,1)
scatter(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),60,'filled')
[cc,pp]=corr(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),'rows','complete','type','Spearman')
title(['Place Field Start: corr=' num2str(cc) 'pvalue='  num2str(pp) ])
xlabel('Theta Scores')
ylabel('Place Field Limit Shift (cm)')

lsline

subplot(1,2,2)
scatter(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),60,'filled')
[cc,pp]=corr(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),'rows','complete','type','Spearman')
title(['Place Field End: corr=' num2str(cc) 'pvalue='  num2str(pp)])
xlabel('Theta Scores')
ylabel('Place Field Limit Shift (cm)')


lsline
linkaxes

%  To DO:
% asse x = posizione del picco 
% asse y: 1) asimmetry
% 2) shift start 
% 3) shift end

 figure(660), clf
 %%%%% PF START and END shift as a function of PF peak
 subplot(1,2,1)
scatter(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),60,'filled')
[cc,pp]=corr(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),'rows','complete','type','Spearman')
title(['PF start shift as a function of PF peak; corr=' num2str(cc) 'pvalue='  num2str(pp)])
lsline
xlabel('PF peak on the track')
ylabel('PF start shift')

 subplot(1,2,2)
scatter(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),60,'filled')
[cc,pp]=corr(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),'rows','complete','type','Spearman')
title(['PF end shift as a function of PF peak; corr=' num2str(cc) 'pvalue='  num2str(pp)])
lsline
xlabel('PF peak on the track')
ylabel('PF end shift')
linkaxes


 %%%%% Asymmetry as a function of PF peak
 sel = find(dspl<remap_value);
 figure(665)
 sgtitle('Asymmetry as a function of PF peak')

 subplot(1,2,1)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak Position')
title('Laser OFF')

subplot(1,2,2)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak Position')
title('Laser ON')
lsline
linkaxes

% Divide precessing e locking 
for cond = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

if(cond==1) %locking
    sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cond==2) %precessing
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end


 figure(6665)
 sgtitle('Asymmetry as a function of PF peak')
 if cond == 1
 subplot(2,2,cond)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Locking Laser OFF')

 subplot(2,2,cond+1)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Locking Laser ON')


 elseif cond == 2
    subplot(2,2,cond+1)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Precessing Laser OFF')

 subplot(2,2,cond+2)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline  
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Precessing Laser ON')
 end 

lsline
linkaxes


end 




%% Cell properties (asymmetry, shifts) plots - WITH VELOCITY CONTROL
close all
clear

%Velocity interval
for Vel =1:4;


direction_name = {'UP','DOWN'};
laser_name = {'OFF','ON'};

cell_name = {'LOCKING','PRECESSING'};

% data_origin = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\!LaserTrondheim\Matteo_GLMLaser_speed\';
 data_origin = 'F:\!papers\!LaserTrack_Paper\!DATA_LaserTrondheim\CorrectData\';
% prop_dir = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\';
% glm_dir = 'C:\Users\Matteo\Desktop\LaserTrondheim\Matteo_GLMLaser\';

gro = 0;

%Chose which cells to take depending of remap of PEAK 
remap_value = 20;

    
theta_All_OFF = [];
theta_All_ON = [];
pf_l_all = [];
pfsasy_All = [];

for Laser = [0 1]

Sparsity_All = [];
Skaggs_All = [];
pfs_All = [];
MeanR_All = [];
MaxR_All = [];
MPos_All = [];




for anm = [1 2 3 4 5 6]
for sess = 1:3
            
   if anm == 1 && sess == 2 || anm == 1 && sess == 3 || anm == 2 && sess == 2 || anm == 2 && sess == 3 || anm == 3 && sess == 3 || anm == 4 && sess == 1 || anm == 5 && sess == 1
   continue   
   end 
    
for direction = 1:2

    
     disp(['Doing Something_' direction_name{direction} '_' laser_name{Laser+1} '_A' num2str(anm) '_S' num2str(sess) '_F' num2str(1)])
    
load([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) 'Sp2Vel' num2str(Vel)],...
    'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')    

% load(['/Users/federico/GitHub/Basins_GLM_paper/Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],...
%     'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')   
    
Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 

if(direction == 2)
   pf_asym = -pf_asym; 
   pf_limits = 100 - fliplr(pf_limits);
   
   
end

pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>0 & pf_limits(:,2)<100)); 

if Laser ==0 
    theta_All_OFF = cat(1,theta_All_OFF,theta_scores);
elseif Laser == 1
    theta_All_ON = cat(1,theta_All_ON,theta_scores);
end 

pf_l_all = cat(1,pf_l_all,pf_limits);

size(theta_scores)
size(MaxRate)


end
end
end 

end


n_cells = size(pf_l_all,1)/2;
pf_1 = pf_l_all(1:n_cells,:);
pf_2 = pf_l_all(n_cells+1:end,:);

pa_1 = pfsasy_All(1:n_cells);
pa_2 = pfsasy_All(n_cells+1:end);


for cond = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

if(cond==1)
sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cond==2)
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end

figure(301)
subplot(4,2,(Vel-1)*2+cond)
histogram(pa_1(sel),-10:0.5:10,'Normalization','probability')
hold on
histogram(pa_2(sel),-10:0.5:10,'Normalization','probability')
xlabel('Asymmetry')
ylabel('Probability')
title(cell_name{cond})
linkaxes
hold on



figure(302)
subplot(4,2,(Vel-1)*2+cond)
[F,X]=ecdf(pa_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pa_2(sel));

plot(X,F,'LineWidth',2)
xlabel('Asymmetry')
ylabel('Cumulative Probability')
h(cond) = kstest2(pa_1(sel),pa_2(sel))
title([cell_name{cond} ': h=' num2str(h(cond))])
linkaxes
hold on

figure(5002)
subplot(4,2,(Vel-1)*2+cond)
[F,X]=ecdf(pf_2(sel,1)-pf_1(sel,1));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pf_2(sel,2)-pf_1(sel,2));
plot(X,F,'LineWidth',2)

[F,X]=ecdf(pf_2(sel,3)-pf_1(sel,3));
plot(X,F,'LineWidth',2)
xlim([-50 50])
plot([0 0],[0 1],'--k')
plot([-50 50],[0.5 0.5],'--k')

xlabel('Shift (cm)')
ylabel('Cumulative Probability')
title(cell_name{cond})
legend({'Start','Centre','End'},'FontSize',16,'Location','northwest')
hold on

end

figure(301)
legend({laser_name{1}, laser_name{2}},'FontSize',16)

figure(302)
legend({laser_name{1}, laser_name{2}},'FontSize',16,'Location','east')



figure(5001)
histogram(pf_2(:,1)-pf_1(:,1),-50:4:50,'Normalization','probability')
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field Start')

histogram(pf_2(:,2)-pf_1(:,2),-50:4:50,'Normalization','probability')
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field Centre')

histogram(pf_2(:,3)-pf_1(:,3),-50:4:50,'Normalization','probability')
xlabel('Shift (cm)')
ylabel('Probability')
title('Place Field End')

linkaxes


figure(6001)
 subplot(4,2,(Vel-1)*2+1)
scatter(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),60,'filled')
xlabel('Theta Scores')
ylabel('Place Field Limit Shift (cm)')
title('Place Field Start')
lsline
hold on
subplot(4,2,(Vel-1)*2+2)
scatter(theta_All_OFF(dspl<remap_value),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),60,'filled')
xlabel('Theta Scores')
ylabel('Place Field Limit Shift (cm)')
title('Place Field End')
lsline
linkaxes
hold on
%%%%%%%%%%%%%%%%%

 figure(660)
 %%%%% PF START and END shift as a function of PF peak
 subplot(4,2,(Vel-1)*2+1)
scatter(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,1)-pf_1(dspl<remap_value,1),60,'filled')
lsline
xlabel('PF peak on the track')
ylabel('PF start shift')
title('PF start shift as a function of PF peak')
hold on
 subplot(4,2,(Vel-1)*2+1)
scatter(pf_2(dspl<remap_value,2),pf_2(dspl<remap_value,3)-pf_1(dspl<remap_value,3),60,'filled')
lsline
xlabel('PF peak on the track')
ylabel('PF end shift')
title('PF end shift as a function of PF peak')
linkaxes
hold on

 %%%%% Asymmetry as a function of PF peak
 sel = find(dspl<remap_value);
 figure(665)
 sgtitle('Asymmetry as a function of PF peak')

 subplot(1,2,1)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak Position')
title('Laser OFF')

subplot(1,2,2)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak Position')
title('Laser ON')
lsline
linkaxes

% Divide precessing e locking 
for cond = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

if(cond==1) %locking
    sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cond==2) %precessing
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end


 figure(6665)
 sgtitle('Asymmetry as a function of PF peak')
 if cond == 1
 subplot(2,2,cond)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Locking Laser OFF')

 subplot(2,2,cond+1)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Locking Laser ON')


 elseif cond == 2
    subplot(2,2,cond+1)
scatter(pf_1(sel,2),pa_1(sel),60,'filled')
lsline
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Precessing Laser OFF')

 subplot(2,2,cond+2)
scatter(pf_2(sel,2),pa_2(sel),60,'filled')
lsline  
ylabel('Asymmetry')
xlabel('PF Peak')
title('Phase Precessing Laser ON')
 end 

lsline
linkaxes


end 
end
