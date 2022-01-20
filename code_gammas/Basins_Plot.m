
clear all
close all

%% Load Data

lay_contact = [1 2 5]; % Position of CA1 layers with respect to probe contacts 


pyr = lay_contact(1)
rad = lay_contact(2)
slm = lay_contact(3)

% Select your data
animal = 1;
session = 1; 


% Load the LFP and CSD signals
load(['Animal' num2str(animal) '\'  num2str(session) 'C_Raw_CSD'],'Raw_CSD')
load(['Animal' num2str(animal) '\'  num2str(session) 'C_Raw_LFP'],'Raw_LFP')

% If you want to compute the CSD from the LFP
clear Raw_CSD
con3 = 1;
for ch=2:size(Raw_LFP,1)-1
   con3 = con3 + 1;
   Raw_CSD(con3-1,:) = - Raw_LFP(con3-1,:)+2*Raw_LFP(con3,:)-Raw_LFP(con3+1,:);
   
end



%% Wavelet Comodulogram 

%Number of channels in your CSD  your recordings!
NChannels = 14;

kk=0;
% Variable for phase amplitudse coupling computed with respect of different
% reference oscillations (e.g. Theta taken from different layers)
Phase_Energy_Lay={};

%Select from which channel to take your reference oscillation (e.g. Theta oscillation)
for pick_th= [slm]

    
kk=kk+1;
% Recording sample rate
srate = 1000;

% Wavelet decomposition for all channels
cwt_all=[];
for ch_i = 1:NChannels
lfp = Raw_CSD(ch_i,:);

% Pick your frequency intervals
Freq_Lim = [2 300];

[cwt_lay,Frq]=cwt(lfp,'amor',srate,'FrequencyLimits',Freq_Lim);

cwt_all=cat(3,cwt_all,cwt_lay);


end

% Frequency interval for the reference oscillation
Carrier_Frq = [6 10];

th_frq=find(Frq<Carrier_Frq(2)  & Frq>Carrier_Frq(1));

Th_Sig=squeeze(sum(cwt_all(th_frq,:,pick_th),1));
amp = abs(Th_Sig);
amp = zscore(amp);

% Phase of reference oscillation
ang = angle(Th_Sig);


% Compute phase amplitude coupling 
Phase_Energy=[];
for ch_i=1:NChannels
All_Amp=squeeze(abs(cwt_all(:,:,ch_i)));
nbin = 36; % number of phase bins
phasebins = -pi:2*pi/nbin:pi;

Ph_Cou=zeros(size(All_Amp,1),nbin);

cc=0;
for ph=phasebins(1:end-1)
cc=cc+1;

tt=find(ang<ph+2*pi/nbin & ang>ph & amp>0);

    Ph_Cou(:,cc)=mean(All_Amp(:,tt),2);
    
end

Phase_Energy=cat(3,Phase_Energy,Ph_Cou);


end
Phase_Energy_Lay{kk}=Phase_Energy;

end


%%Basins - Construct 

nbin = 36; % number of phase bins
phasebins = -pi:2*pi/nbin:pi;
% Just repeat the phase cycle 3 times for visualization 
window_bin = [phasebins(1:end-1)-2*pi phasebins(1:end-1) phasebins(1:end-1)+2*pi];


comod_str=zeros(3,3);
comod_pha=zeros(3,3);
comod_lay=zeros(3,3);
comod_frq=zeros(3,3);

% Cycle through different reference oscillations as computed above
for Layer = 1:1 

Phase_Energy=Phase_Energy_Lay{Layer};

% Rearrange frequency range
Frq_U=flipud(Frq);
Phase_Energy=Phase_Energy(end:-1:1,:,:);

frq_t = find(Frq_U<200 & Frq_U>15);
Frq_U=Frq_U(frq_t);
Phase_Energy=Phase_Energy(frq_t,:,:);

Ener_Peaks_All=zeros(size(Phase_Energy,1)*size(Phase_Energy,2)*3,14);
norm = Phase_Energy(:,:,ch_i);
norm = min(norm,[],2);
norm = repmat(norm,1,size(Phase_Energy,2)*3);

EE_Cut_All=[];
for ch_i = 1:NChannels


Energy=repmat(Phase_Energy(:,:,ch_i),1,3);
%EE_N=Energy./norm;
EE=(imgaussfilt(Energy./norm,[1,1]));
EE_N=(imgaussfilt(Energy./norm,[1.5,1.5]));
EE_Z = zscore(EE_N,[],2);


% Compute Laplacian 
LLV=zeros(size(EE_N));
LLO=zeros(size(EE_N));
for cc=1:size(EE_N,2)
LLV(:,cc)=del2(EE_N(:,cc));
end

for rr=1:size(EE_N,1)
LLO(rr,:)=del2(EE_N(rr,:));
end

% Take phase-frequency regions with negative curvature 
LLO(LLO>=0)=0;
LLO(LLO<0)=1;
LLV(LLV>=0)=0;
LLV(LLV<0)=1;


LL_Tot = LLO.*LLV;

LL_Tot(LL_Tot<=0)=0;
LL_Tot(LL_Tot>0)=1;

% Select connected negative-curvature regions
Freq_Field = bwconncomp(LL_Tot);


EE_Cut=zeros(size(EE));

% Find significantly coupled bins in negtive-curvature regions
for ff=1:Freq_Field.NumObjects
pix = Freq_Field.PixelIdxList{ff};
    pz=find(EE_Z(pix)>1);


    EE_Cut(pix(pz)) = EE_N(pix(pz));
end
EE_Cut_All=cat(3,EE_Cut_All,EE_Cut);


figure(1)
subplot(131)
imagesc(EE_N)
title('Phase-Amplitude Coupling')
xlabel('Frequency')
ylabel('Phase')
set(gca,'YDir','reverse')
subplot(132)
imagesc(LL_Tot)
title('Negative Curvature Regions')
xlabel('Frequency')
ylabel('Phase')
set(gca,'YDir','reverse')
subplot(133)
imagesc(EE_Cut)
title('Significant Regions')
xlabel('Frequency')
ylabel('Phase')
set(gca,'YDir','reverse')
pause(0.1)


end

% Total 3D basin organization across layers and frequencies 
Basins_Struct{Layer} = EE_Cut_All;


EE_Cut_All = Basins_Struct{Layer}; 
v = EE_Cut_All(:,:,:);
% v=EE_Cut_All(:,nbin+1:2*nbin,:);
min_x = 0;
max_x = 10;
min_y = 0;
max_y = 10;

figure(2001)

for ch=1:14

    imgzposition = ch;
planeimg = v(:,:,ch)';
    surf([min_x max_x],[min_y max_y],repmat(imgzposition, [2 2]),...
    planeimg,'facecolor','texture','FaceAlpha','texture',...
        'AlphaDataMapping','direct',...
        'AlphaData',10+planeimg*80)

hold on

end
caxis([0 1.5])
%alpha 0.5
set(gca,'YDir','reverse')
set(gca,'ZDir','reverse')






%% Basins - CoModulation Maximum and Position

% Extract gamma specific coupling strength and phase

v=EE_Cut_All(:,nbin+1:2*nbin,:);


% Limits for gamma ranges

G_Lim = [20 45; 60 90; 100 180];

for gamma = 1:3

Frq_G=find(Frq_U>G_Lim(gamma,1) & Frq_U<G_Lim(gamma,2));

v_G = squeeze(max(max(v(Frq_G,:,1:9),[],3),[],1));

[comod_str(Layer, gamma), comod_pha(Layer, gamma)] = max(v_G);

v_G = squeeze(max(max(v(Frq_G,:,1:9),[],2),[],1));

[comod_str(Layer, gamma), comod_lay(Layer, gamma)] = max(v_G);

v_G = squeeze(max(max(v(Frq_G,:,1:9),[],2),[],3));

[~,cc] = max(v_G);
comod_frq(Layer, gamma) = Frq_U(Frq_G(cc));

end






%% Basins - 3D plotting 

region_size_limit = 40;

DC=distinguishable_colors(100);
figure(50 + Layer)
clf;

y=1:size(v,1);
x=phasebins(1:end-1);
z=1:size(v,3);
ll=0;
for lev=[0.999 1.01 1.05]
    ll=ll+1;
p = patch(isosurface(x,y,z,v,lev));
hold on
isonormals(x,y,z,v,p)
p.FaceColor = DC(ll,:);
p.EdgeColor = 'none';
alpha(0.5)
daspect([1 1 1])
view(3); 
axis square

end

% Take 3D connected significant bins
Field_Chann=bwconncomp(EE_Cut_All,6);

figure(80 + Layer)
clf;
cc=0;

% Cycle through connected regions
for comp=1:Field_Chann.NumObjects
    
    % Take only regions comprising more than x bins 
    if(numel(Field_Chann.PixelIdxList{comp})>region_size_limit)
    cc=cc+1;
        v=zeros(size(EE_Cut_All));
    v(Field_Chann.PixelIdxList{comp})=1;%EE_Cut_All(Field_Chann.PixelIdxList{comp});
    
            c_v=zeros(size(EE_Cut_All));
    c_v(Field_Chann.PixelIdxList{comp})=EE_Cut_All(Field_Chann.PixelIdxList{comp});
    
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

ylim([20  250])
daspect([1 1 1])
view(3); 
axis square
camlight 
lighting gouraud
grid on 


%% Basins - Single Gamma Plotting

figure(90 + Layer)
clf;

for gamma = 1:3
subplot(1,3,gamma)

Frq_G=find(Frq_U>G_Lim(gamma,1) & Frq_U<G_Lim(gamma,2));
v_mask = zeros(size(EE_Cut_All));
v_mask(Frq_G,:,:)=1;



for comp=1:Field_Chann.NumObjects
    
    if(numel(Field_Chann.PixelIdxList{comp})>region_size_limit)
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

ylim([20  250])
daspect([1 1 1])
view(3); 
axis square
camlight 
lighting gouraud
grid on 








end





end
