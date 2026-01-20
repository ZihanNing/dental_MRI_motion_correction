
function [stateSample] = stateSampleFromRec(rec)

k = cat(1, rec.Assign.z{2}, rec.Assign.z{3});
k = k.';

NY = size(rec.y);

%%% Shift to make compatible with DISORDER recon
kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
%k(:,1) = (NY(2) - (k(:,1)-1)) - kShift(2); % 2nd PE 
%k(:,2) = (NY(3) - (k(:,2)-1)) - kShift(3); % 3rd PE = slices

k(:,1) = -(k(:,1) + kShift(2) - NY(2) ) + 1 ; % 2nd PE 
k(:,2) = -(k(:,2) + kShift(3) - NY(3) ) + 1 ; % 2nd PE 


%%% Flip sampling and offset (because Siemens fft vs. ifft convention)
% if mod(NY(2),2)==0;offset =-1;else;offset=0;end
% Lines = mod( Lines-1 + offset, NY(2) ) +1;
% if mod(NY(3),2)==0;offset =-1;else;offset=0;end
%Partitions = mod( Partitions -1 + offset, NY(3) ) +1;

if mod(NY(2),2)==0;offset =-1;else;offset=0;end
k(:,1) = mod(k(:,1)-offset -1, NY(2)) +1 ;
if mod(NY(3),2)==0;offset =-1;else;offset=0;end
k(:,2) = mod(k(:,2) - offset-1, NY(3))  +1;

%% Detect linear direction
etl = zeros([1 2]);
th = .1;
for i=1:2
    a = k(:,i) ;

    shotStart =  find(abs(diff(a))> round(th*max(a)) );

    if isempty(shotStart)
        stateSample = 1;
    else
        stateSample = makeBinsStart(a, shotStart);
    end
    
    figure; 
    subplot 311; plot(a); hold on ; scatter(shotStart, a(shotStart))
	subplot 312; plot(diff(a)); hold on ; scatter(shotStart, dynInd(diff(a),shotStart,1))
    subplot 313; plot(stateSample);
    
    etl(i) =   mean( accumarray(stateSample(:),onesL(stateSample)) );
end

%     kTrajs=abs(fct(k-mean(k,1))); %YB: By doing this frequency analysis, you don't actually need the sampling patter procided to the scanner!
%     N=size(kTrajs);
%     [~,iAs]=max(dynInd(kTrajs,1:floor(N(1)/2),1),[],1);    
%     iM=min(iAs,[],2);
%     NSamples=size(k,1)/(iM-1);
%     stsa = makeBins(size(k,1), iM-1);
    
[~,dim] = min( abs(etl - rec.Par.Labels.RFEchoTrainLength));

%% Final one
a = k(:,dim) ;
shotStart =  find(abs(diff(a))> round(th*max(a)) );
stateSample = makeBinsStart(a, shotStart);

    