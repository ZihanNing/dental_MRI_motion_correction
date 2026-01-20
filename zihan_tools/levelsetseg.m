function Mask_Vol = levelsetseg(Vol, varargin)
%LEVELSETSEG Slice-wise DRLSE segmentation for 3D volumes
%
% Mask_Vol = levelsetseg(Vol)
% Mask_Vol = levelsetseg(Vol, 'preset', 'MPRAGE')
% Mask_Vol = levelsetseg(Vol, 'preset', 'PDwSPACE')
% Mask_Vol = levelsetseg(Vol, 'lambda',..., 'alfa',...)
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 09-Jan-2026

%% -------------------- defaults --------------------
param.lambda   = 3.0;
param.alfa     = 1.5;
param.epsilon  = 3;
param.sigma    = 2;

param.timestep   = 5;
param.mu         = 0.2 / param.timestep;
param.iter_inner = 5;
param.iter_outer = 100;
param.initial_msk_rows = 5:60;
param.initial_msk_cols = 5:60;

preset = '';

%% -------------------- parse varargin --------------------
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'preset'
            preset = varargin{k+1};
        case 'lambda'
            param.lambda = varargin{k+1};
        case 'alfa'
            param.alfa = varargin{k+1};
        case 'epsilon'
            param.epsilon = varargin{k+1};
        case 'sigma'
            param.sigma = varargin{k+1};
        case 'timestep'
            param.timestep = varargin{k+1};
            param.mu = 0.2 / param.timestep;
        case 'iter_inner'
            param.iter_inner = varargin{k+1};
        case 'iter_outer'
            param.iter_outer = varargin{k+1};
        otherwise
            error('Unknown parameter: %s', varargin{k});
    end
end

%% -------------------- presets --------------------
if ~isempty(preset)
    switch lower(preset)
        case 'mprage'
            param.lambda  = 3.0;
            param.alfa    = 1.5;
            param.epsilon = 3;
            param.sigma   = 2;
            param.initial_msk_rows = 5:60;
            param.initial_msk_cols = 5:60;

        case 'pdwspace'
            param.lambda  = 3.5;
            param.alfa    = 4;
            param.epsilon = 3;
            param.sigma   = 1.2;
            param.initial_msk_rows = 3:61;
            param.initial_msk_cols = 3:61;
            
        case 't2wspace'
            param.lambda  = 3.5;
            param.alfa    = 4;
            param.epsilon = 1.5;
            param.sigma   = 1.2;
            param.initial_msk_rows = 3:61;
            param.initial_msk_cols = 3:61;

        otherwise
            error('Unknown preset: %s', preset);
    end
end

%% -------------------- preprocessing --------------------
Vol = double(Vol);
[nx, ny, nz] = size(Vol);
Mask_Vol = false(nx, ny, nz);

rTh=[3 1.2 1.2];
ellip_msk = gather(getEllipsMask(Vol,rTh)==1);

%% -------------------- slice-wise DRLSE --------------------
fprintf('Segmentation by level-set (DRLSE)...... \n')

for i = 1:nz
    
    
    % compute slice by slice
    Img = Vol(:,:,i);
    Cal_Col = max(size(Vol,1),64);Cal_Row = max(size(Vol,2),64);
    Img = imresize(Img, [Cal_Col Cal_Row], 'bilinear');

    timestep=5;  % time step
    mu=0.2/timestep;  % coefficient of the distance regularization term R(phi),0.2
    iter_inner=5;
    iter_outer=100;
    lambda=param.lambda; % coefficient of the weighted length term L(phi),5
    alfa=param.alfa;  % coefficient of the weighted area term A(phi), 1.5
    epsilon=param.epsilon; % papramater that specifies the width of the DiracDelta function, 1.5

    sigma=param.sigma;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.

    % initialize LSF as binary step function
    c0=2;
    initialLSF=c0*ones(size(Img));
    % generate the initial region R0 as a rectangle
    initialLSF(param.initial_msk_rows,param.initial_msk_cols)=-c0;  
    phi=initialLSF;
%     initialLSF = c0 - 2*c0 * ellip_msk;
%     phi=initialLSF;
    

%     figure(1);
%     mesh(-phi);   % for a better view, the LSF is displayed upside down
%     hold on;  contour(phi, [0,0], 'r','LineWidth',2);
%     title('Initial level set function');
%     view([-80 35]);

%     figure(2);
%     imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%     title('Initial zero level contour');
%     pause(0.5);

    potential=2;  
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
end
