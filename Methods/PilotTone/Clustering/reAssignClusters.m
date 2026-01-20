
function [stateSampleNew, XClusteredNew, dClustersNew, dSamplesNew] = reAssignClusters(X, stateSample, XClustered, dClusters, dSamples)

%REASSIGNCLUSTERS re-assigns samples to a set of clusters so that each cluster approximately has the same number of samples.
%
%   [STATESAMPLENEW,XCLUSTERED,DCLUSTERS,DSAMPLES] = REASSIGNCLUSTERS(X,STATESAMPLE,XCLUSTERED,DCLUSTERS,DSAMPLES)
%
%   Yannick Brackenier 2023-07-10
%
%TODO: change clusters centers when you start re-arranging : Need to do a
%linear weighting with where samples comes from

%%% SET PARAMETERS
fillBiggestFirst = 1;
NClust = max(stateSample);
NSamples = length(stateSample);

%%% INITIALISE
newBinIdx = makeBins(NSamples, NClust);
nSa = accumarray(newBinIdx(:), onesL(newBinIdx(:)));%Samples to assign to each cluster
if fillBiggestFirst; [~,orderClust] = sort(nSa,'descend');else;[~,orderClust] = sort(nSa,'acend');end
assert(multDimSum(nSa)==NSamples,'reAssignClusters:: Total number of samples does not add up.');

stateSampleNew = inf*onesL(stateSample);

%%% GO OVER CLUSTERS AND RE-ASSIGN
distClust = computeDistance(X,XClustered);%[NSamples x NClust]
for k=1:NClust
    %Extract cluster
    clustId = orderClust(k);
    temp = distClust(:,clustId);
    %Sort
    [~,sampleIdx] = sort(temp(:),'ascend');
    sampleIdx = sampleIdx(1:nSa(k));
    stateSampleNew(sampleIdx) = clustId;
    %Set to inf
    distClust(sampleIdx,:) = inf;
end
assert(~any(isinf(stateSampleNew)),'reAssignClusters:: Total number of samples does not add up.');

%%% COMPUTE NEW CLUSTER
dClustersNew = inf*onesL(dClusters);
XClusteredNew = inf*onesL(XClustered);
for k=1:NClust
    idx = stateSampleNew==k;
    XClusteredNew(k,:) = multDimMea(X(idx,:),1);
    dClustersNew(k) = multDimMea(X(idx,:)-XClusteredNew(k,:),1:2); 
end
dSamplesNew = computeDistance(X,XClusteredNew);

end


%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%
function [dSamp] = computeDistance(X,XClustered)
    %X has size [NSamples x NFeatures]
    
    %Get the number of clusters
    K = size(XClustered,1);
    
    %Compute the distance per sample
    dSamp = zeros(size(X,1),K,'like',X); 
    for k=1:K
        dSamp(:,k)= multDimMea(abs(X(:,:) - XClustered(k,:)),2);
    end
end




