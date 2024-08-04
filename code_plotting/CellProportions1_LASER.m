
%%Plotting and calculating PF differences - NO VELOCITY CONTROL

close all

clear

direction_name = {'UP','DOWN'};
laser_name = {'OFF','ON'};
cell_name = {'LOCKING','PRECESSING'};

%data_origin = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\!LaserTrondheim\Matteo_GLMLaser_speed\';
 data_origin = 'F:\!papers\!LaserTrack_Paper\!DATA_LaserTrondheim\CorrectData\';
 %data_origin = 'F:\processing\!Laser_HD\control_animals\CTRL_Data_2\'

% prop_dir = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\';
% glm_dir = 'C:\Users\Matteo\Desktop\LaserTrondheim\Matteo_GLMLaser\';

gro = 0;

%Chose which cells to take depending of remap of PEAK (based on average remapping values)
remap_value = 15;

    
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

Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 

if Laser ==0 
    theta_All_OFF = cat(1,theta_All_OFF,theta_scores);
elseif Laser == 1
    theta_All_ON = cat(1,theta_All_ON,theta_scores);
end 

if(direction == 2)
   pf_asym = -pf_asym; 
   pf_limits = 100 - fliplr(pf_limits);
   
   
end
pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>0 & pf_limits(:,2)<100)); 



pf_l_all = cat(1,pf_l_all,pf_limits);

size(theta_scores)
size(MaxRate)


end
end
end

if Laser == 0
    Skaggs_All_OFF = Skaggs_All 
    MeanR_All_OFF = MeanR_All
    MaxR_All_OFF = MaxR_All
    pf_size_OFF = pfs_All
    MPos_OFF = MPos_All
elseif Laser == 1
    Skaggs_All_ON = Skaggs_All
    MeanR_All_ON = MeanR_All
    MaxR_All_ON = MaxR_All
    pf_size_ON = pfs_All
    MPos_ON = MPos_All
end

end


n_cells = size(pf_l_all,1)/2;
pf_1 = pf_l_all(1:n_cells,:);
pf_2 = pf_l_all(n_cells+1:end,:);
%
pa_1 = pfsasy_All(1:n_cells);
pa_2 = pfsasy_All(n_cells+1:end);
%
Skaggs_1 = Skaggs_All_OFF
Skaggs_2 = Skaggs_All_ON
%
MeanR_1 = MeanR_All_OFF;
MeanR_2 = MeanR_All_ON;
%
MaxR_1 = MaxR_All_OFF;
MaxR_2 = MaxR_All_ON;
%
pf_size_1 = pf_size_OFF
pf_size_2 = pf_size_ON
%
baricenter_1 = MPos_OFF;
baricenter_2 = MPos_ON;


%test pf center

for cell_type =1:2
dspl=abs(pf_1(:,2)-pf_2(:,2)); % based on peak remapping
% dspl = pf_1(:,2); % based on pf center

total_remap = (find(dspl>remap_value))
no_remap = (find(dspl<remap_value))


if(cell_type==1)
%  sel = find(theta_All_OFF<0 & dspl<89 & dspl>11);% based on pf center
    sel = find(theta_All_OFF<0 & dspl<remap_value);

elseif(cell_type==2)
%  sel = find(theta_All_OFF>0 & dspl<89 & dspl>11);% based on pf center
     sel = find(theta_All_OFF>0 & dspl<remap_value);
end

% figure(10)
% subplot(1,2,cell_type)
% [F,X]=ecdf(Skaggs_1(sel));
% plot(X,F,'LineWidth',2)
% hold on
% [F,X]=ecdf(Skaggs_2(sel));
% plot(X,F,'LineWidth',2)
% 
% xlabel('Skaggs Info')
% ylabel('Cumulative Probability')
% [h(cell_type),p(cell_type)] = ttest(Skaggs_1(sel),Skaggs_2(sel))
%  title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
% linkaxes
% legend ('Laser OFF','Laser ON')
% sgtitle('Skaggs information')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(20)
subplot(1,2,cell_type)
[F,X]=ecdf(MeanR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MeanR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Mean rate')
ylabel('Cumulative Probability')
[h(cell_type),p(cell_type)] = ttest(MeanR_1(sel),MeanR_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Mean FR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(30)
subplot(1,2,cell_type)
[F,X]=ecdf(MaxR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MaxR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Max rate')
ylabel('Cumulative Probability')
[h(cell_type),p1(cell_type)] = ttest(MaxR_1(sel),MaxR_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p1(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Max FR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(40)
subplot(1,2,cell_type)
[F,X]=ecdf(pf_size_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pf_size_2(sel));
plot(X,F,'LineWidth',2)

xlabel('pf size')
ylabel('Cumulative Probability')
[h(cell_type),p2(cell_type)] = ttest(pf_size_1(sel),pf_size_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p2(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Place field size')

end 

figure(1500)
subplot(1,3,1)
scatter(theta_All_OFF,theta_All_ON,60,'filled')
lsline
xlabel('Theta Score OFF')
ylabel('Theta Score ON')
[cc,pp]=corr(theta_All_OFF(:),theta_All_OFF(:),'rows','complete','type','Spearman')
title(['Theta Score OFF vs ON, all cells: pval=' num2str(cc)])
refline(1,0)
subplot(1,3,2)
scatter(theta_All_OFF(total_remap),theta_All_ON(total_remap),60,'filled')
lsline
xlabel('Theta Score OFF')
ylabel('Theta Score ON')
[cc,pp]=corr(theta_All_OFF(total_remap),theta_All_OFF(total_remap),'rows','complete','type','Spearman')
title(['Theta Score OFF vs ON, REMAP cells: pval=' num2str(cc)])
refline(1,0)
subplot(1,3,3)
scatter(theta_All_OFF(no_remap),theta_All_ON(no_remap),60,'filled')
lsline
xlabel('Theta Score OFF')
ylabel('Theta Score ON')
[cc,pp]=corr(theta_All_OFF(no_remap),theta_All_OFF(no_remap),'rows','complete','type','Spearman')
title(['Theta Score OFF vs ON, NO REMAP cells: pval=' num2str(cc)])
refline(1,0)
% 
% figure(1550)
% sgtitle('Theta Score Distribution OFF vs ON')
% histogram(theta_All_OFF,-0.4:0.1:0.6,'Normalization','probability')
% hold on
% histogram(theta_All_ON,-0.4:0.1:0.6,'Normalization','probability')
% legend('Theta Score OFF','Theta Score ON')
% linkaxes

figure(50)
histogram(baricenter_1,0:2:40,'Normalization','probability')
hold on 
histogram(baricenter_2,0:2:40,'Normalization','probability')
xlabel('Baricenter')
ylabel('Probability')
bc = mean(baricenter_1)
bc2 = mean(baricenter_2)
title(['Baricenter OFF: ' num2str(bc) '; Baricenter ON: ' num2str(bc2) ])

figure(60)
histogram(dspl,0:5:90,'Normalization','probability')
xlabel('PF Peak shift in cm')
ylabel('Probability')
md = mean(dspl)
title(['Remapping mean value =' num2str(md) 'cm'])


%% Plotting and calculating PF differences - WITH VELOCITY CONTROL

close all

clear

for Vel = 1:4

direction_name = {'UP','DOWN'};
laser_name = {'OFF','ON'};
cell_name = {'LOCKING','PRECESSING'};

%  data_origin = 'F:\!papers\!LaserTrack_Paper\!DATA_LaserTrondheim\CorrectData\';
 data_origin = 'F:\processing\!Laser_HD\control_animals\Cell_PropsAllNew\';

% data_origin = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\!LaserTrondheim\Matteo_GLMLaser_speed\';
% prop_dir = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\';
% glm_dir = 'C:\Users\Matteo\Desktop\LaserTrondheim\Matteo_GLMLaser\';

gro = 0;

% %Velocity interval
% Vel = 1;

%Chose which cells to take depending of remap of PEAK 
remap_value = 22;

    
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

%     
%     if(~exist([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) '.mat'],'file'))
%     continue
%     end
    
     disp(['Doing Something_' direction_name{direction} '_' laser_name{Laser+1} '_A' num2str(anm) '_S' num2str(sess) '_F' num2str(1)])
    
load([data_origin 'c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) 'Sp2Vel' num2str(Vel)],...
    'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')    

% 
% load([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) 'Sp2Vel' num2str(Vel)],...
%    'theta_scores')    

       
% load(['/Users/federico/GitHub/Basins_GLM_paper/Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro)],...
%     'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')   
    
Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 

if Laser ==0 
    theta_All_OFF = cat(1,theta_All_OFF,theta_scores);
elseif Laser == 1
    theta_All_ON = cat(1,theta_All_ON,theta_scores);
end 

if(direction == 2)
   pf_asym = -pf_asym; 
   pf_limits = 100 - fliplr(pf_limits);
   
   
end
pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>0 & pf_limits(:,2)<100)); 



pf_l_all = cat(1,pf_l_all,pf_limits);

size(theta_scores)
size(MaxRate)


end
end
end

if Laser == 0
    Skaggs_All_OFF = Skaggs_All 
    MeanR_All_OFF = MeanR_All
    MaxR_All_OFF = MaxR_All
    pf_size_OFF = pfs_All
    
elseif Laser == 1
    Skaggs_All_ON = Skaggs_All
    MeanR_All_ON = MeanR_All
    MaxR_All_ON = MaxR_All
    pf_size_ON = pfs_All
    
end

end


n_cells = size(pf_l_all,1)/2;
pf_1 = pf_l_all(1:n_cells,:);
pf_2 = pf_l_all(n_cells+1:end,:);
%
pa_1 = pfsasy_All(1:n_cells);
pa_2 = pfsasy_All(n_cells+1:end);
%
Skaggs_1 = Skaggs_All_OFF
Skaggs_2 = Skaggs_All_ON
%
MeanR_1 = MeanR_All_OFF;
MeanR_2 = MeanR_All_ON;
%
MaxR_1 = MaxR_All_OFF;
MaxR_2 = MaxR_All_ON;
%
pf_size_1 = pf_size_OFF
pf_size_2 = pf_size_ON
%


clear sel 
for cell_type = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

if(cell_type==1)
    sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cell_type==2)
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end

% figure(10)
% % subplot(1,2,cell_type)
% subplot(4,2,(Vel-1)*2+cell_type)
% [F,X]=ecdf(Skaggs_1(sel));
% plot(X,F,'LineWidth',2)
% hold on
% [F,X]=ecdf(Skaggs_2(sel));
% plot(X,F,'LineWidth',2)
% 
% xlabel('Skaggs Info')
% ylabel('Cumulative Probability')
% clear h
% [h(cell_type),p(cell_type)] = ttest(Skaggs_1(sel),Skaggs_2(sel))
%  title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
% linkaxes
% legend ('Laser OFF','Laser ON')
% sgtitle('Skaggs information')
% hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(20)
subplot(4,2,(Vel-1)*2+cell_type)
[F,X]=ecdf(MeanR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MeanR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Mean rate')
ylabel('Cumulative Probability')
clear h
[h(cell_type),p(cell_type)] = ttest(MeanR_1(sel),MeanR_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Mean FR')
hold on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(30)
subplot(4,2,(Vel-1)*2+cell_type)
[F,X]=ecdf(MaxR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MaxR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Max rate')
ylabel('Cumulative Probability')
clear h
[h(cell_type),p(cell_type)]= ttest(MaxR_1(sel),MaxR_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Max FR')
hold on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(40)
subplot(4,2,(Vel-1)*2+cell_type)
[F,X]=ecdf(pf_size_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pf_size_2(sel));
plot(X,F,'LineWidth',2)

xlabel('pf size')
ylabel('Cumulative Probability')
clear h
[h(cell_type),p(cell_type)]= ttest(pf_size_1(sel),pf_size_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend ('Laser OFF','Laser ON')
sgtitle('Place field size')
hold on

end 
end 

figure(66)
histogram(dspl,0:5:90,'Normalization','probability')
xlabel('PF Peak shift in cm')
ylabel('Probability')
md = mean(dspl)
title(['Remapping mean value =' num2str(md) 'cm'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control analysis with Alt1 e Alt2

close all
clc
clear

direction_name = {'UP','DOWN'};
laser_name = {'OFF','ON'};
cell_name = {'LOCKING','PRECESSING'};

% data_origin = 'C:\Users\Matteo\Desktop\LaserTrack_Paper\!LaserTrondheim\Matteo_GLMLaser_speed\';
 data_origin = 'F:\!papers\!LaserTrack_Paper\!DATA_LaserTrondheim\CorrectData\';

%data_origin = '~/Dropbox/Projects_NIJ/Matteo/Trondheim_LinearLaser/';
gro =0;
pf_limits_all={};

gro = 0;

% %Velocity interval
% Vel = 1;

%Chose which cells to take depending of remap of PEAK 
remap_value = 15;

    
theta_All_OFF = [];
theta_All_ON = [];
pf_l_all = [];
pfsasy_All = [];

%FISSO IL LASER E PLOTTO DIFFERENZE ALT1 e ALT2
for alt = [1 2]
for Laser =1 

Sparsity_All = [];
Skaggs_All = [];
pfs_All = [];
MeanR_All = [];
MaxR_All = [];
MPos_All = [];
pfsasy_All = [];
theta_All = [];
pf_limits_temp=[];

for anm = [1 2 3 4 5 6]
for sess = 1:3

    
for direction = 1:2

    
   if anm == 1 && sess == 2 || anm == 1 && sess == 3 || anm == 2 && sess == 2 || anm == 2 && sess == 3 || anm == 3 && sess == 3 || anm == 4 && sess == 1 || anm == 5 && sess == 1
   continue   
   end   
    
     disp(['Doing Something_' direction_name{direction} '_' laser_name{Laser+1} '_A' num2str(anm) '_S' num2str(sess) '_F' num2str(1)])
    
load([data_origin 'Cell_Props/c_props_A' num2str(anm) 'S' num2str(sess) 'D' num2str(direction) 'G' num2str(gro) 'L' num2str(Laser) 'Sp2Alt' num2str(alt)],...
    'Sparsity','Skaggs','MeanRate','MaxRate','pf_size','Max_Pos','pf_limits','pf_asym','theta_scores')    
    

if (direction == 2)
pf_asym = -pf_asym
pf_limits = 100 - fliplr(pf_limits);
end


Sparsity_All= cat(1,Sparsity_All,Sparsity);
Skaggs_All = cat(1,Skaggs_All,Skaggs);
pfs_All= cat(1,pfs_All,pf_size);
MeanR_All = cat(1,MeanR_All,MeanRate);
MaxR_All = cat(1,MaxR_All,MaxRate);
MPos_All= cat(1,MPos_All,Max_Pos); 
theta_All = cat(1,theta_All,theta_scores);

if alt == 1
theta_All_OFF = cat(1,theta_All_OFF,theta_scores);
end 

size(theta_scores)
size(MaxRate)

pfsasy_All = cat(1,pfsasy_All,pf_asym(pf_limits(:,2)>0 & pf_limits(:,2)<100)); 
% pf_l_all = cat(1,pf_l_all,pf_limits);
pf_limits_temp = cat(1,pf_limits_temp,pf_limits);

end
end
end 


HH = histcounts(MPos_All,0:40,'Normalization','probability');
HH = smoothdata(HH,'movmean',5);


figure(10)
histogram(pfsasy_All,-10:1:10,'Normalization','probability')
hold on
PF_asy(Laser+1)= mean(abs(pfsasy_All))
title('PF Asymemtry')
figure(30)
histfit(Skaggs_All,35,'normal')
title('Skaggs Info')
hold on 
figure(40)
histfit(MaxR_All,35,'normal')
title('Max rate')
hold on 
figure(50)
histfit(MeanR_All,35,'normal')
title('Mean rate')
hold on 
figure(1)
histogram(pfs_All)
title('PF Size')
hold on


PF_size(alt)= mean(pfs_All)
Max_FR(alt)= mean(MaxR_All)
Mean_FR(alt)= mean(MeanR_All)
Sparsity_Info(alt)= mean(Sparsity_All)
Skaggs_Info(alt)= mean(Skaggs_All)
Baricenter(:,alt) = MPos_All;

% if alt == 1
%     pf_limits_all{alt} =pf_limits_temp ;
% elseif alt == 2
%     pf_limits_all{alt} =pf_limits_temp ;
% end 

end

if alt == 1
    Skaggs_All_alt1 = Skaggs_All 
    MeanR_All_alt1 = MeanR_All
    MaxR_All_alt1 = MaxR_All
    pf_size_alt1 = pfs_All
    pf_limits_temp1 = pf_limits_temp
    Baricenter1 = Baricenter
    
elseif alt == 2
    Skaggs_All_alt2 = Skaggs_All
    MeanR_All_alt2 = MeanR_All
    MaxR_All_alt2 = MaxR_All
    pf_size_alt2 = pfs_All  
    pf_limits_temp2 = pf_limits_temp
    Baricenter2 = Baricenter
end

end


n_cells = size(pf_limits_temp1,1);
 pf_1 = pf_limits_temp1
 pf_2 = pf_limits_temp2
%
pa_1 = pfsasy_All(1:n_cells);
pa_2 = pfsasy_All(n_cells+1:end);
%
Skaggs_1 = Skaggs_All_alt1
Skaggs_2 = Skaggs_All_alt2
%
MeanR_1 = MeanR_All_alt1;
MeanR_2 = MeanR_All_alt2;
%
MaxR_1 = MaxR_All_alt1;
MaxR_2 = MaxR_All_alt2;
%
pf_size_1 = pf_size_alt1
pf_size_2 = pf_size_alt2
%


clear sel 

for cell_type = 1:2
dspl=abs(pf_1(:,2)-pf_2(:,2));

total_remap = numel(find(dspl>remap_value))
no_remap = numel(find(dspl<remap_value))

if(cell_type==1)
    sel = find(theta_All_OFF<0 & dspl<remap_value);
elseif(cell_type==2)
    sel = find(theta_All_OFF>0 & dspl<remap_value);
end

% figure(10)
% subplot(1,2,cell_type)
% [F,X]=ecdf(Skaggs_1(sel));
% plot(X,F,'LineWidth',2)
% hold on
% [F,X]=ecdf(Skaggs_2(sel));
% plot(X,F,'LineWidth',2)
% 
% 
% xlabel('Skaggs Info')
% ylabel('Cumulative Probability')
% h(cell_type) = kstest2(Skaggs_1(sel),Skaggs_2(sel))
%  title([cell_name{cell_type} ': h=' num2str(h(cell_type))])
% linkaxes
% legend('Odd laps', 'Even laps')
% sgtitle('Skaggs information')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(20)
subplot(1,2,cell_type)
[F,X]=ecdf(MeanR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MeanR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Mean rate')
ylabel('Cumulative Probability')
[h(cell_type), p(cell_type)] = ttest(MeanR_1(sel),MeanR_2(sel))
title([cell_name{cell_type} ': h=' num2str(p(cell_type)) '; h=' num2str(h(cell_type))])
linkaxes
legend('Odd laps', 'Even laps')
sgtitle('Mean FR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(30)
subplot(1,2,cell_type)
[F,X]=ecdf(MaxR_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(MaxR_2(sel));
plot(X,F,'LineWidth',2)

xlabel('Max rate')
ylabel('Cumulative Probability')
[h(cell_type), p(cell_type)]  = ttest(MaxR_1(sel),MaxR_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend('Odd laps', 'Even laps')
sgtitle('Max FR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(40)
subplot(1,2,cell_type)
[F,X]=ecdf(pf_size_1(sel));
plot(X,F,'LineWidth',2)
hold on
[F,X]=ecdf(pf_size_2(sel));
plot(X,F,'LineWidth',2)

xlabel('pf size')
ylabel('Cumulative Probability')
[h(cell_type), p(cell_type)]  = ttest(pf_size_1(sel),pf_size_2(sel))
 title([cell_name{cell_type} ': h=' num2str(p(cell_type))])
linkaxes
legend('Odd laps', 'Even laps')
sgtitle('Place field size')



end 
 

figure(50),clf
histogram(dspl,0:5:90,'Normalization','probability')
xlabel('PF Peak shift in cm')
ylabel('Probability')
md = mean(dspl)
title(['Remapping mean value =' num2str(md) 'cm'])

figure(60),clf
histogram(Baricenter1,0:2:40,'Normalization','probability')
hold on 
histogram(Baricenter2,0:2:40,'Normalization','probability')
xlabel('Baricenter')
ylabel('Probability')
bc = mean(Baricenter1)
bc2 = mean(Baricenter2)
title(['Baricenter OFF: ' num2str(bc) '; Baricenter ON: ' num2str(bc2) ])
legend ('Laser OFF','Laser ON')