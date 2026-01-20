function ACS2rec(fileAcq,writeRAWFlag,writeNIIFlag)

%%% LOAD THE REC OF HIGH-RES SCAN
if ~iscell(fileAcq)%Single acquisition
    [pathOuTemp,file]=fileparts(fileAcq);
    if exist(strcat(fileAcq,'.mat'),'file')~=2; dat2Rec(fileAcq);end
    ss= load(strcat(fileAcq,'.mat'));rec=ss.rec;ss=[]; 
else %Multiple acquisitions
    disp('Error: not supported for multiple acquisitions currently!')
    pause;
end

%%% GENERATE recS STRUCTURE 
recS = rec;
recS.y = rec.ACS;
recS.Enc.AcqVoxelSize = rec.Enc.UnderSampling.ACSVoxelSize;
recS.Par.Mine.APhiRec = rec.Par.Mine.APhiACS;

recS.Names.pathOu = pathOuTemp;
recS.Names.Name = file;
recS.Plan.Suff=''; recS.Plan.SuffOu='';

%%% ZN: FILTER
filterRefData = 0
if filterRefData
    gibbsRing=.5;
    %Create filter
    sp = 1;
    NS=multDimSize(recS.y,1:3);
    HS=buildFilter(2*NS,'tukeyIso',sp,0,gibbsRing,1);
    %Filter
    recS.y=filtering(recS.y,HS,1);     
end

%%% COIL SENSITIVITY MAP ESTIMATION
doDivEstimation=0;
if doDivEstimation
%     recS=solveSensit7T_Div(recS, []); % ZN: currently not support recB, recB = []
    recS=solveSensit7T_Div_noresignresol(recS, []); % ZN: currently not support recB, recB = []
else
    tic
%     recS=solveSensit7T(recS, []); % ZN: currently not support recB, recB = []
%     recS=solveSensit7T_lucilio(recS, []); % ZN: debugging
    recS=solveSensit7T_noresignresol(recS, []); % ZN: debugging
    disp('Time for espirit:(in minutes)')
    espirit_time = toc/60
end

%%% SAVE RAW DATA
recS.doDivEstimation = doDivEstimation;
recS.ACSflag = 1; % ZN: use ACS line for estimation
if writeRAWFlag
    if existsFileVar(strcat(pathOuTemp,file),'rec')==2
        recS.y=[];recS.N=[];recS.Assign=[];%Remove raw data to avoid storing twice
    end
    save(strcat(pathOuTemp,file,'_ACS.mat'),'recS');%If this file is also used for something else, don't overwrite .mat file
    fprintf('Raw file saved (recS by ACS line):\n   %s\n', strcat(file,'_ACS.mat'));
end

%%% SAVE NIFTI
if writeNIIFlag>0
    fprintf('Writing NIFTI files.\n');
    pathOuNII = strcat(pathOuTemp, '/Re-Se/'); if ~exist( pathOuNII,'dir');mkdir(pathOuNII);end
    fileSave = strcat( pathOuNII,file);
    xW=[];xW{1} = recS.S; 
    MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; 
    MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
    writeNII(fileSave, {'ACS_Se'},xW, MSW, MTW);
end
if writeNIIFlag>0
    fprintf('Writing NIFTI files.\n');
    pathOuNII = strcat(pathOuTemp, '/Re-Se/'); if ~exist( pathOuNII,'dir');mkdir(pathOuNII);end
    fileSave = strcat( pathOuNII,file);
    xW=[];xW{1} = recS.x; 
    MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; 
    MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
    writeNII(fileSave, {'ACS_rec'},xW, MSW, MTW);
end
tmp = RSOS(recS.S); % ZN: to see if there's any odd estimation 
if writeNIIFlag>0
    fprintf('Writing NIFTI files.\n');
    pathOuNII = strcat(pathOuTemp, '/Re-Se/'); if ~exist( pathOuNII,'dir');mkdir(pathOuNII);end
    fileSave = strcat( pathOuNII,file);
    xW=[];xW{1} = tmp; 
    MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; 
    MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
    writeNII(fileSave, {'ACS_coilmap'},xW, MSW, MTW);
end

end

