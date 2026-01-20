
function [label, elementsFound] = inRange (x, rr, typeEdge, typeInt)

%INRANGE finds the elements of an array that are within a set of numerical intervals and labels each element accordingly.
%   [IDX,ELEMENTSFOUND]=INRANGE(X, RR, {TYPEEDGE},{TYPEINT})
%   * X is the array to look at.
%   * RR is a list with the intervals defining the ranges to look for. The number of intervals equal length(rr)-1.
%   * {TYPEEDGE} determines how to use the edges of the extreme intervals  when labeling the elements in X.
%             'includeboth' : Both end points included in interval (default)
%             'includeleft' : Left end point only included in interval
%   * {TYPEINT} determines how to use the edges of the non-extreme (read: middle) intervals  when labeling the elements in X.
%             'includeright': Right end point only included in interval
%             'excludeboth' : Neither end point included in interval
%   ** LABEL is an array of size size(X) with each element having the index of interval it belongs to. Elements not belonging to any interval have index 0.
%   ** ELEMENTSFOUND is the number of elements that belong to any interval.
%
%   Yannick Brackenier 2022-08-30

if nargin<2 || isempty(rr); [~,~,rr] = getRange(x);end
if nargin<3 || isempty(typeEdge); typeEdge = 'includeboth';end 
if nargin<4 || isempty(typeInt); typeInt = 'includeleft';end %either 'includeleft' or 'includeright'

NInt = length(rr)-1;
label = zerosL(real(x));
elementsFound = zeros([1 NInt]);

for i=1:NInt
    bound = rr(i:(i+1));
    
    if i==1 %%% Start interval
        if strcmp(typeEdge,'includeboth') || strcmp(typeEdge,'includeleft')
            idxTemp = bound(1)<=x(:);
        else
            idxTemp = bound(1)<x(:);
        end
        
        if strcmp(typeInt,'includeleft') && ~(NInt==1 && (strcmp(typeEdge,'includeboth') || strcmp(typeEdge,'includeright')) )%make sure that itworks for NInt==1 as well
            idxTemp = idxTemp & (x(:)<bound(2));
        else 
            idxTemp = idxTemp & (x(:)<=bound(2));
        end
    elseif i==NInt %%% End interval
        
        if strcmp(typeInt,'includeleft')
            idxTemp = x(:)>=bound(1);
        else 
            idxTemp = x(:)>bound(1);
        end
        
        if strcmp(typeEdge,'includeboth') || strcmp(typeEdge,'includeright')
            idxTemp = idxTemp & (x(:)<=bound(2));
        else
            idxTemp = idxTemp & (x(:)<bound(2));    
        end
        
    else %%% Middle intervals
        if strcmp(typeInt,'includeleft')
            idxTemp = (bound(1)<= x(:)) & (x(:)<bound(2));
        elseif strcmp(typeInt,'includeright')
            idxTemp = (bound(1)< x(:)) & (x(:)<=bound(2));
        end
    end
    label(idxTemp) = i;
    elementsFound(i) = multDimSum(idxTemp);
end

%elementsFound = any(label(:));
