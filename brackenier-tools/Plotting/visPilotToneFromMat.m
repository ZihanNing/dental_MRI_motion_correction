function []= visPilotToneFromMat(fileName)

load(strcat(fileName,'.mat'));

p = PT.pTimeTest;
if size(p,3)>1; p=dynInd(p,PT.idxMB,3);end
visPTSignal(p,[],[],[], 1, [], [], [], 'Original');

visPTSignal(p./sqrt(normm(p,[],1)),[],[],[], 1, [], [], [], 'Normalised');