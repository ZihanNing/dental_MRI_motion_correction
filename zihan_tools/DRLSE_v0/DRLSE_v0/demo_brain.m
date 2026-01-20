%  This Matlab code demonstrates an edge-based active contour model as an application of 
%  the Distance Regularized Level Set Evolution (DRLSE) formulation in the following paper:
%
%  C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", 
%     IEEE Trans. Image Processing, vol. 19 (12), pp. 3243-3254, 2010.
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/

clear all;
close all;

load('/home/gadgetron/matlab/3D-moco-siemens/zihan_tools/meas_MID00031_FID143573_PDwSPACE_orig_sepACS_FAvar_DISORDER_move_ACS.mat');
Vol=double(abs(recS.x));
Vol = permute(Vol,[2 3 1]);
figure,imshow3D(Vol,[0 2000])
Num_Slice = size(Vol,3);
%% parameter setting
% Img = Img(:,:,32);
% Img = imresize(Img, [64 64], 'bilinear');
Mask_Vol = zeros(size(Vol));
tic
for i = 1:Num_Slice
    
    % compute slice by slice
    Img = Vol(:,:,i);
    Img = imresize(Img, [64 64], 'bilinear');

    timestep=5;  % time step
    mu=0.2/timestep;  % coefficient of the distance regularization term R(phi),0.2
    iter_inner=5;
    iter_outer=100;
    lambda=3.0; % coefficient of the weighted length term L(phi),5
    alfa=1.5;  % coefficient of the weighted area term A(phi), 1.5
    epsilon=4; % papramater that specifies the width of the DiracDelta function, 1.5

    sigma=2;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.

    % initialize LSF as binary step function
    c0=2;
    initialLSF=c0*ones(size(Img));
    % generate the initial region R0 as a rectangle
    initialLSF(5:60, 5:60)=-c0;  
    phi=initialLSF;

%     figure(1);
%     mesh(-phi);   % for a better view, the LSF is displayed upside down
%     hold on;  contour(phi, [0,0], 'r','LineWidth',2);
%     title('Initial level set function');
%     view([-80 35]);

%     figure(2);
%     imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%     title('Initial zero level contour');
%     pause(0.5);

    potential=3;  
    if potential ==1
        potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
    elseif potential == 2
        potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
    else
        potentialFunction = 'double-well';  % default choice of potential function
    end


    % start level set evolution
    for n=1:iter_outer
        phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
        if mod(n,2)==0
%             figure(2);
%             imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
        end
    end

    % refine the zero level contour by further level set evolution with alfa=0
    alfa=0;
    iter_refine = 10;
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);

%     finalLSF=phi;
%     figure(2);
%     imagesc(Img,[0, 2000]); axis off; axis equal; colormap(gray); 
%     hold on;  contour(phi, [0,0], 'r'); 
%     % phi_thres = max(phi(:)) - 0.2*(max(phi(:)) - min(phi(:)));
%     % hold on; contour(phi, [phi_thres,phi_thres], 'g');
%     str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
%     title(str);

%     pause(1);
%     figure;
%     mesh(-finalLSF); % for a better view, the LSF is displayed upside down
%     hold on;  contour(phi, [0,0], 'r','LineWidth',2);
%     str=['Final level set function, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
%     title(str);
%     axis on;
    
    % Get the mask
    mask = phi<=0;
    Mask_Vol(:,:,i) = imresize(mask,[size(Vol,1),size(Vol,2)],'nearest');
    
    fprintf('Completed %d/%d. \n',i,size(Vol,3))

end
toc

%%
figure,imshow3D(abs(Vol.*Mask_Vol),[0 2000])
figure,imshow3D(abs(Mask_Vol),[])

%%
figure,imshow3D(abs(Vol.*Mask_Vol_polished),[0 2000])
figure,imshow3D(abs(Mask_Vol_polished),[])

