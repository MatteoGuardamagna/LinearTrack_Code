
%% Gamma bands and coupling with Theta
close all
%Routines
sr = 1000;
Freq_Lim = [2 300];
[~,Frq]=cwt(ones(10000,1),'amor',sr,'FrequencyLimits',Freq_Lim);
Frq_U=flipud(Frq);
Frq_U=Frq_U(22:end-5);

%Number of Channels
NChannels = 14;
nbin = 36;

%Layers of your animals
lay_n = [3 6 9;
    2 5 8;
    1 6 8;
    1 5 8;
    2 5 9;
    1 6 9];


for gamma = [1 2 3]

gg=0;
kk=0;

co_str_all = NaN(3,10,6,3);

for animal = [1 2 3 4 5 6]
    
kk=kk+1;
comod_frq_all=[]; 
comod_lay_all=[];
comod_str_all=[];
comod_pha_all=[];
co_str = NaN(3,10,3);

ss_i = 0;
for ss=1:10
    ss_i = ss_i+1;
    if animal == 1 && ss<1
        continue
    elseif animal == 5 && ss>9
        continue
    else
        
   %Loading your data
   load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinData_A' num2str(animal) '_S' num2str(ss) '.mat'])
    end
    
for Layer_Th = 1:3   
   EE_Cut_All = Basins_Struct{Layer_Th}; 
   
   v=EE_Cut_All(:,nbin+1:2*nbin,:);


G_Lim = [20 45; 60 90; 100 180];


Frq_G=find(Frq_U>G_Lim(gamma,1) & Frq_U<G_Lim(gamma,2));

v_G = squeeze(max(max(v(Frq_G,:,:),[],1),[],2));

co_str(:,ss_i,Layer_Th) = v_G(lay_n(animal,:));


end
   
%    comod_frq_all = cat(3,comod_frq_all,comod_frq);
%    comod_lay_all = cat(3,comod_lay_all,comod_lay);
%    comod_str_all = cat(3,comod_str_all,comod_str);
%    comod_pha_all = cat(3,comod_pha_all,comod_pha);
    
end

for Layer_Th = 1:3
figure(90+gamma)
gg = gg+1;
subplot(6,3,gg)

dd = smoothdata(co_str(:,:,Layer_Th),2,'movmean',2);

co_str_all(:,:,animal,Layer_Th) = dd; 


plot(dd','LineWidth',3)
%ylim([1 1.7])
ylim([1 2])
hold on 
%plot(max(dd,[],1),'--k')

Cou_Str(animal,Layer_Th,gamma,:) = (mean(dd(:,:),2));
Cou_Str_z(animal,Layer_Th,gamma,:) = zscore(mean(dd(:,:),2));
end
% 
% for ff = 1:3
%     gg= gg+1;
%     subplot(4,3,gg)
%     
% plot(squeeze(comod_str_all([1 2 3],ff,:))')
% end

end


figure(100)
for ll = 1:3
    subplot(3,3,(gamma-1)*3+ll)
errorbar(nanmean(co_str_all(:,:,:,ll),3)',(nanstd(co_str_all(:,:,:,ll),[],3)')/2.5,'LineWidth',3)
%ylim([1 1.7])
ylim([1 2])
hold on 
end

end

figure(100)
subplot(3,3,1)
title('Pyr Theta')
ylabel('Slow Gamma')
subplot(3,3,2)
title('Rad Theta')
subplot(3,3,3)
title('SLM Theta')
subplot(3,3,4)
ylabel('Medium Gamma')
subplot(3,3,7)
ylabel('Fast Gamma')

legend

sgtitle({'Gamma Bands (rows, from the 3 Layers, b=pyr; r=rad; y=slm)','Coupling with Theta (columns, from the 3 Layers)'})

%%
close all

lay_n = [3 6 9;
    2 5 8;
    1 6 8;
    1 5 8;
    2 5 9;
    1 6 9];
Layer = 2;
gg=0;
kk=0;

%For Animal
for animal = [1 2 3 4 5 6]
kk=kk+1;
comod_frq_all=[];
comod_lay_all=[];
comod_str_all=[];
comod_pha_all=[];

for ss=4:10
   load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinData_A' num2str(animal) '_S' num2str(ss) '.mat'])
   comod_frq_all = cat(3,comod_frq_all,comod_frq);
   comod_lay_all = cat(3,comod_lay_all,comod_lay);
   comod_str_all = cat(3,comod_str_all,comod_str);
   comod_pha_all = cat(3,comod_pha_all,comod_pha);

end

for ff = 1:3
    gg= gg+1;
    subplot(4,3,gg)
    
plot(squeeze(comod_lay_all([1 2 3],ff,:))')
end


end

%%
close all

for typ = 1:2

% Define your Layers
lay_n = [3 6 9;
    2 5 8;
    1 6 8;
    1 5 8;
    2 5 9;
    1 6 9];

ph_mea = zeros(3,10,7);
ph_off = zeros(1,7);

va_mea = zeros(3,10,7);

kk=0;
for animal = [1 2 3 4 5]
    
sr = 1000;
Freq_Lim = [2 300];
[~,Frq]=cwt(ones(10000,1),'amor',sr,'FrequencyLimits',Freq_Lim);
Frq_U=flipud(Frq);
Frq_U=Frq_U(22:end-5);

NChannels = 14;
nbin = 36; % number of phase bins
phasebins = -pi:2*pi/nbin:pi;
window_bin = [phasebins(1:end-1)-2*pi phasebins(1:end-1) phasebins(1:end-1)+2*pi];

if(animal<6)
s_list = [1:2];
else
    s_list = [1:1];
end
for ss=s_list
kk=kk+1;
l_pl = 0;
for gamma = 1:3
    l_pl = l_pl+1;
%for Layer = 1:3

G_Lim = [20 45; 60 90; 100 180];
%gamma = 2;
Layer = 3; 
ph_all=0;

if(typ==1)
    
   load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinsData_RUN_A' num2str(animal) '_S' num2str(ss) '.mat'])
elseif(typ==2)
    load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinData_REM_A' num2str(animal) '_S' num2str(ss) '.mat'])

end


EE_Cut_All = Basins_Struct{Layer}; 

v=EE_Cut_All(:,nbin+1:2*nbin,:);
%% Basins - 3D plotting 
DC=distinguishable_colors(100);

Frq_G=find(Frq_U>G_Lim(gamma,1) & Frq_U<G_Lim(gamma,2));
v_mask = zeros(size(EE_Cut_All));
v_mask(Frq_G,:,:)=1;

[v_i,ph_i]=max(squeeze(max(v.*(v_mask(:,nbin+1:2*nbin,:)),[],1)),[],1);

ph_i = unwrap(phasebins(ph_i));

%ph_i = ph_i(lay_n(animal,:));

ph_i = interp1(1:14,ph_i,(linspace(lay_n(animal,1),lay_n(animal,3),7)));
v_i = interp1(1:14,v_i,(linspace(lay_n(animal,1),lay_n(animal,3),7)));
ph_i = smoothdata(ph_i,'movmean',1);
if(gamma==1)
    ph_ref = ph_i;
elseif(gamma==2)
    ph_diff = ph_i-ph_ref;
    ph_diff = exp(1i*ph_diff);
     
end

ph_i = exp(1i*ph_i);


ph_all = ph_all + ph_i;




figure(70)
subplot(1,3,l_pl)
plot((angle(ph_all)),1:7)
hold on
set(gca,'YDir','reverse')
%xlim([-3.14 3.14])


ph_mea(l_pl,kk,:) = ph_all;%ph_mea(l_pl,:) + ph_all;
va_mea(l_pl,kk,:) = v_i;%va_mea(l_pl,:) + v_i;
if(gamma==2)
ph_off = ph_off + ph_diff;
end
end

end
end
%ph_mea(2,1)=ph_mea(2,2)+0.2;
figure(71)
subplot(1,2,typ)
errorbar(squeeze(circ_mean(angle(ph_mea),[],2))',squeeze(circ_std(angle(ph_mea),[],[],2))'/3,'LineWidth',2)
xticks([1,4,7])
xticklabels({'Pyr','Rad','SLM'})
ylabel('Theta Phase')
legend('Slow','Medium' ,'Fast')
ylim([-pi , pi])
if (typ == 1)
    title('RUN')
elseif(typ == 2)
title('REM')
end


figure(72)
subplot(1,2,typ)
plot(unwrap(angle(ph_off)))

figure(73)
subplot(1,2,typ)


errorbar(squeeze(mean(va_mea,2))',squeeze(std(va_mea,[],2))'./3,'LineWidth',2)

xticks([1,4,7])
xticklabels({'Pyr','Rad','SLM'})
ylabel('CouplingStrength')
legend('Slow','Medium' ,'Fast')
ylim([0 2])


if (typ == 1)
    title('RUN')
elseif(typ == 2)
title('REM')





end
end

%% Plot Basins 

close all

for typ = [1,2]


kk=0;
for animal = [1]
    kk=kk+1;
sr = 1000;
Freq_Lim = [2 300];
[~,Frq]=cwt(ones(10000,1),'amor',sr,'FrequencyLimits',Freq_Lim);
Frq_U=flipud(Frq);
Frq_U=Frq_U(22:end-5);

NChannels = 14;
nbin = 36; % number of phase bins
phasebins = -pi:2*pi/nbin:pi;
window_bin = [phasebins(1:end-1)-2*pi phasebins(1:end-1) phasebins(1:end-1)+2*pi];

Layer = 3;

G_Lim = [20 45; 60 90; 100 180];
%gamma = 1;

ss =1; 
for gamma=[1 2 3]
    if typ==1
    load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinsData_RUN_A' num2str(animal) '_S' num2str(ss) '.mat'])
    elseif typ==2
    load(['./OldGLM_Basins/AnimalNG' num2str(animal) '/BasinData_REM_A' num2str(animal) '_S' num2str(ss) '.mat'])
    end


EE_Cut_All = Basins_Struct{Layer}; 

v=EE_Cut_All(:,nbin+1:2*nbin,:);
%% Basins - 3D plotting 
DC=distinguishable_colors(100);

Frq_G=find(Frq_U>G_Lim(gamma,1) & Frq_U<G_Lim(gamma,2));
v_mask = zeros(size(EE_Cut_All));
v_mask(Frq_G,:,:)=1;



Field_Chann=bwconncomp(EE_Cut_All,6);

figure(80 + Layer)
subplot(1,2,typ)
%figure()
%clf;
cc=0;
for comp=1:Field_Chann.NumObjects
    
    if(numel(Field_Chann.PixelIdxList{comp})>40)
    cc=cc+1;
        v=zeros(size(EE_Cut_All));
    v(Field_Chann.PixelIdxList{comp})=1;%EE_Cut_All(Field_Chann.PixelIdxList{comp});
    
    v = v.*v_mask;
    
            c_v=zeros(size(EE_Cut_All));
    c_v(Field_Chann.PixelIdxList{comp})=EE_Cut_All(Field_Chann.PixelIdxList{comp});
    
    c_v = c_v.*v_mask;
    
    Fr_Ch_Val=max(reshape(c_v,[],NChannels),[],1);
    
y=Frq_U;
x=window_bin;
z=1:size(v,3);
p = patch(isosurface(x,y,z,v,0.99));
isonormals(x,y,z,v,p)


Col_chann=zeros(size(p.Vertices,1),1);

for ch_c=1:NChannels

    Col_chann(find(round(p.Vertices(:,3))==ch_c))=Fr_Ch_Val(ch_c);

end

p.FaceVertexCData=Col_chann;

p.FaceColor = 'interp';
p.EdgeColor = 'none';

    end

end

zlim([0 16])
set(gca,'ztick',[1:14])
set(gca, 'ZDir','reverse')

ylim([23  200])
xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
daspect([1 1 1])
view([49 7]); 
axis square
camlight 
lighting gouraud
grid on 
caxis([1 2])

%pause()
end


end
end

