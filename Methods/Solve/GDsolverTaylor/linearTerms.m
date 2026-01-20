
function [K1, K2, K3] = linearTerms(N, scaling, dimStore)
    
    if nargin < 2 || isempty(scaling); scaling = ones(size(N));end
    if nargin < 3 || isempty(dimStore); dimStore=4;end
    assert(dimStore>3, 'linearTerms:: dimStore should be bigger than 3 to avoid conflict with spatial dimensions.');
    
    NCent=centerIdx(N);
    
    %%% create linear terms  
    rangex = (1:N(1)) -  NCent(1); rangex = rangex / (N(1)/2); %ranges from -1 -> 1
    rangey = (1:N(2)) -  NCent(2); rangey = rangey / (N(2)/2);
    rangez = (1:N(3)) -  NCent(3); rangez = rangez / (N(3)/2);
    
    rangex = scaling(1) * rangex;%ranges from -scale -> scale
    rangey = scaling(2) * rangey;
    rangez = scaling(3) * rangez;
    
    test = ifftshift(rangex); 
    assert(test(1)==0,'linearTerms:: Zero-frequency term of the linear term not correct.');
    
    [Kx, Ky, Kz] = ndgrid(rangex, rangey, rangez);
    
    %%% return
    if nargout==1
        K1 = cat(dimStore, Kx, Ky, Kz);
    elseif nargout==3
        K1=Kx; K2=Ky;K3=Kz;
    else
        error('linearTerms:: No valid output/number of outputs') 
    end
        
end