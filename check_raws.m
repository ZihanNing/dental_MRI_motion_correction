clc;
clear;
close all;

%% Add path
addpath(genpath('/home/gadgetron/matlab/Reconstructionv4/mapVBVD'))
addpath(genpath('/home/gadgetron/matlab/usual_used/'))
addpath(genpath('/home/gadgetron/matlab/3D-moco-siemens'))

%% Read raw
TW = dat2TWIX('meas_MID00627_FID09434_FLAIR_TwinsUK.dat');
TW = TW{2}; % do not care the noise data
TW.image.flagIgnoreSeg=1;%See mapVBVD tutorial
TW.image.dataSize(1:11) % RO, CHA, PE, 3D

% get the kspace to be showed
data =dynInd(TW.image,{1:TW.image.dataSize(2) ':'},[2 length(TW.image.dataSize)]); 
size(data) %RO-Channel-PE1-PE2
data = squeeze(data(:,30,:,:,1,1)); % pick one channel and one average/contrast data

% permute the data to proper view [3D, PE]
data = permute(data,[3 2 1]); % 3D, PE, RO
figure,imshow3D(abs(data),[0 2e-5])

%% show the image
figure;
for i = 1:16
    slab_to_show = 100+19*i;
    subplot(4,4,i), imshow(abs(data(:,:,slab_to_show)),[0 2e-5]), title(['RO direction ',num2str(slab_to_show)]);
end
%% show the signal profile of one specific region
data = permute(data,[3,2,1]); % RO, Lin, Par
% % DISORDER
% region_early = 122:124; % Lin
% region_later = 126:128; % Lin
% region_centre = 100:102; % Lin

% linear
region_early = 126:128; % Lin
region_later = 1:3; % Lin
region_centre = 100:102; % Lin

% profile in the RO direction
RO_axis = 1:TW.image.dataSize(1);
RO_pro = zeros(3,length(RO_axis));
for i = RO_axis
    profile = data(i,region_early,:); RO_
    pro(1,i) = abs(sum(profile,'all'))./numel(profile); % early echoes
    profile = data(i,region_centre,:); RO_pro(2,i) = abs(sum(profile,'all'))./numel(profile); % centre
    profile = data(i,region_later,:); RO_pro(3,i) = abs(sum(profile,'all'))./numel(profile); % later echoes
end
figure;
subplot(2,2,1),plot(RO_axis,RO_pro(1,:),'r');
hold on; plot(RO_axis,RO_pro(3,:),'b');
title('RO profile'),xlabel('RO'),ylabel('averaged signal intensity [mag]'),legend('early echo','later echo')
subplot(2,2,2),plot(RO_axis,RO_pro(2,:));
title('RO profile'),xlabel('RO'),ylabel('averaged signal intensity [mag]'),legend('centre kspace')

% profile in the Par direction
Par_axis = 1:TW.image.dataSize(4);
Par_pro = zeros(3,length(Par_axis));
for i = Par_axis
    profile = data(:,region_early,i); Par_pro(1,i) = abs(sum(profile,'all'))./numel(profile); % early echoes
    profile = data(:,region_centre,i); Par_pro(2,i) = abs(sum(profile,'all'))./numel(profile); % centre
    profile = data(:,region_later,i); Par_pro(3,i) = abs(sum(profile,'all'))./numel(profile); % later echoes
end
subplot(2,2,3),plot(Par_axis,Par_pro(1,:),'r');
hold on; plot(Par_axis,Par_pro(3,:),'b');
title('Par profile'),xlabel('Par'),ylabel('averaged signal intensity [mag]'),legend('early echo','later echo')
subplot(2,2,4),plot(Par_axis,Par_pro(2,:));
title('Par profile'),xlabel('Par'),ylabel('averaged signal intensity [mag]'),legend('centre kspace')



