%% Basins Plotting - Common Properties

close all
data_origin = 'I:\Matteo\processing\Norwegian_processing';

for laser = 1:2

%Layers of your animals

%Layers
lay_n = [1 3 6;
 2 6 10
 2 5 9
 3 6 10
 6 9 12 
 1 2 4]

% %Layers
% lay_n = [1 3 6;
%  2 6 10
%  2 5 9
%  3 6 10
%  6 9 12 
%  1 2 4]

ph_mea = zeros(3,10,7);
ph_off = zeros(1,7);

va_mea = zeros(3,10,7);

kk=0;
for animal = [1 2 3 4 5 6]
    
sr = 1000;
Freq_Lim = [2 300];
[~,Frq]=cwt(ones(10000,1),'amor',sr,'FrequencyLimits',Freq_Lim);
Frq_U=flipud(Frq);
Frq_U=Frq_U(22:end-5);

NChannels = 14;
nbin = 36; % number of phase bins
phasebins = -pi:2*pi/nbin:pi;
window_bin = [phasebins(1:end-1)-2*pi phasebins(1:end-1) phasebins(1:end-1)+2*pi];



for ss=1:3
kk=kk+1;
l_pl = 0;
for gamma = 1:2
    l_pl = l_pl+1;
%for Layer = 1:3

G_Lim = [20 45; 60 90; 100 150];

Layer =3; 
ph_all=0;

    if(~exist([ 'I:\Matteo\processing\Norwegian_processing\Animal' num2str(animal) '\' num2str(ss) 'Laser.mat'],'file'))
    continue
    end
        
   %Loading your data
   if laser == 1
   load([ 'BasinDataCSD_OFF_A' num2str(animal) '_S' num2str(ss) '.mat'])
   elseif laser ==2
   load([ 'BasinDataCSD_ON_A' num2str(animal) '_S' num2str(ss) '.mat'])
   end

   
%       %Loading your data
%    if laser == 1
%    load([data_origin '/BasinDataLFP_OFF_A' num2str(animal) '_S' num2str(ss) '.mat'])
%    elseif laser ==2
%    load([data_origin '/BasinDataLFP_ON_A' num2str(animal) '_S' num2str(ss) '.mat'])
%    end
%    
   
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
subplot(1,2,laser)
errorbar(squeeze(circ_mean(angle(ph_mea),[],2))',squeeze(circ_std(angle(ph_mea),[],[],2))'/3,'LineWidth',2)
xticks([1,4,7])
xticklabels({'Pyr','Rad','SLM'})
ylabel('Theta Phase')
legend('Slow','Medium' ,'Fast')
ylim([-pi , pi])
if (laser == 1)
    title('Laser OFF')
elseif(laser == 2)
title('Laser ON')
end


figure(72)
subplot(1,2,laser)
plot(unwrap(angle(ph_off)))

figure(73)
subplot(1,2,laser)
errorbar(squeeze(mean(va_mea,2))',squeeze(std(va_mea,[],2))'./3,'LineWidth',2)

xticks([1,4,7])
xticklabels({'Pyr','Rad','SLM'})
ylabel('CouplingStrength')
legend('Slow','Medium' ,'Fast')
ylim([0 4])
linkaxes

if (laser == 1)
title('Laser OFF')
elseif(laser == 2)
title('Laser ON')


end
end
