

% clear all;
% 
% addpath(genpath('/home/ybr19/Projects/B0-Shimming/Functions'));
% addpath(genpath('/home/ybr19/Projects/Reconstruction'));
% addpath(genpath('/home/ybr19/Projects/PilotTone'));
% addpath(genpath('/home/ybr19/Projects/Simulations/Reconstruction/Methods/synthesizeKSpace'))
% 
% addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'));
% addpath(genpath('/home/ybr19/Software/Utilities'));
% addpath(genpath('/home/ybr19/Software/Registration/rigidReg'));

function  createHybridDMC(studyIdList)

for studyId = studyIdList    
    %% Data path
    if studyId==3%B.E-G.
        studies_20231103_MoPT_Study3_HV3
    elseif studyId==4%D.S.
        studies_20231106_MoPT_Study3_HV4
    elseif studyId==5%Y.X.
        studies_20231123_MoPT_Study3_HV5;
    elseif studyId==6%V.F.
        studies_20231123_MoPT_Study3_HV6;
    elseif studyId==7%J.A.
        studies_20231123_MoPT_Study3_HV7;   
    elseif studyId==8%O.Y.
        studies_20231201_MoPT_Study3_HV8;
    elseif studyId==9%Z.N.
        studies_20231201_MoPT_Study3_HV9;
    elseif studyId==10%P.D.C.
        studies_20231201_MoPT_Study3_HV10;
    elseif studyId==11%M.M.
        studies_20231204_MoPT_Study3_HV11;
    elseif studyId==12%
        studies_20231208_MoPT_Study3_HV12;
    elseif studyId==13%
        studies_20231208_MoPT_Study3_HV13;
    elseif studyId==14%A.T.
        studies_20231208_MoPT_Study3_HV14;
    end
    
    idLowRes = [1:2:11];
    idLowRes = idLowRes(1:6);

    resol = '3mm';

    NAcq = length(idLowRes);
    suff = sprintf('%dscans',NAcq);

    useIntList = [0 1];
    
    %% Load all the low-res ones
    ss = [];
    idx = [];
    for i=1:length(idLowRes)
        fprintf('Loading scan %d\n',i);

        ss{i} = load(strcat(pathIn{1},fileIn{1}{idLowRes(i)},'.mat'));
        ss{i}.yFlat = ss{i}.rec.y;
        ss{i}.pFlat = ss{i}.rec.PT.pSliceImage;
        ss{i}.pTime = ss{i}.rec.PT.pTimeTest;

        %sstemp = load(strcat(pathIn{1},fileIn{1}{idLowRes(i)},'_PT2.mat'));
        %ss{i}.pTime = cat(1, ss{i}.pTime,sstemp.rec.PT.pTimeTest);

        [~, idx{i}] = PESamplesFromRecPTTest(ss{i}.rec);

        % reshape k-space
        N = size( ss{i}.rec.y);
        for n=2:3
            ss{i}.yFlat=fftshiftOperator(ss{i}.yFlat,2,1,n);
            ss{i}.yFlat=fftGPU(ss{i}.yFlat,n)/N(n);
            ss{i}.yFlat=fftshiftGPU(ss{i}.yFlat,n);
        end
        ss{i}.yFlat = resPop(ss{i}.yFlat, 2:3,[],2);

        % reshape PT
        for n=2:3
            ss{i}.pFlat=fftshiftOperator(ss{i}.pFlat,2,1,n);
            ss{i}.pFlat=fftGPU(ss{i}.pFlat,n)/N(n);
            ss{i}.pFlat=fftshiftGPU(ss{i}.pFlat,n);
        end
        ss{i}.pFlat = resPop(ss{i}.pFlat, 2:3,[],2);

    end

    %pSliceImage
    %pTimeTest

    %% Flatten all into vector
    yAll = [];
    idxAll = [];
    kAll =[];
    pTimeAll = [];
    pAll = [];

    for i=1:length(idLowRes)
        yAll = cat(2, yAll, dynInd(ss{i}.yFlat,idx{i},2) );
        pAll = cat(2, pAll, dynInd(ss{i}.pFlat,idx{i},2) );
        kAll = cat(1, kAll, cat( 2, ss{i}.rec.Assign.z{2}(:),  ss{i}.rec.Assign.z{3}(:)));
        offset = 0;%length(idxAll)
        idxAll = cat(1, idxAll, idx{i}(:) + offset);
        pTimeAll = cat(2, pTimeAll, ss{i}.pTime);
    end
    NSamplesTot = size(kAll,1);

    %% Run over perc
    percFullScan = 100/NAcq;
    %assert(mod(percFullScan,1)==0,'');
    percList = [1:NAcq]*percFullScan;
    
    for useInt = useIntList
        for perc = percList
            
            NFiles = perc/percFullScan;
            if NFiles==1 && useInt==1
               shiftList = [0:(NAcq-1)];
               assert(round(perc)==17,'');
            else
               shiftList = 0; 
            end
            
            for shift=shiftList
                %Get samples to extract depending on interl
                if useInt
                    NInterleaves = NAcq;
                    NSamplesPerScan = NSamplesTot/NAcq;
                    idWithinScan = makeBins(NSamplesPerScan, NInterleaves);

                    tt = [];
                    for rep = 0:(NFiles-1)
                        for int=1:NAcq
                           ttTemp = find(   idWithinScan==(mod(int-1+rep+shift,NInterleaves)+1)   ) ;
                           offset = (int-1)*NSamplesPerScan;
                           tt = cat(2, tt, ttTemp+offset); 
                        end
                    end
                    figure ; plot(tt)
                else
                    tt = 1:(NSamplesTot*perc/100);
                    tt = round(tt);
                    offset = shift*NSamplesPerScan;
                    tt = tt + offset;
                    assert(shift==0,'')
                end

                %Extract data of interest
                kTemp = dynInd(kAll, tt,1);
                yTemp = dynInd(yAll, tt,2);
                idxTemp = dynInd(idxAll, tt,1);
                pTemp = dynInd(pAll, tt,2);
                pTimeTemp = dynInd(pTimeAll, tt,2);

                %Split into files
                rec = ss{1}.rec;
                rec.Assign.z{2} = [];
                rec.Assign.z{3} = [];
                rec.y = [];
                rec.PT.pSliceImage = [];
                rec.PT.pTimeTest = [];

                idFile = makeBins(size(yTemp,2), NFiles);
                for i=1:NFiles
                    yTempFile =  dynInd(yTemp,idFile==i ,2); 
                    pTempFile =  dynInd(pTemp,idFile==i ,2); 
                    kTempFile =  dynInd(kTemp,idFile==i ,1);
                    idxTempFile =  dynInd(idxTemp,idFile==i ,1);
                    pTimeTempFile = dynInd(pTimeTemp,idFile==i ,2); 

                    rec.Assign.z{2} = cat(2, rec.Assign.z{2}, resPop(kTempFile(:,1),1,[],2) );
                    rec.Assign.z{3} = cat(2, rec.Assign.z{3}, resPop(kTempFile(:,2),1,[],2) );

                    %Reshape back into single k-space
                    tt = zerosL( resPop(ss{1}.yFlat, 2:3,[],2) );
                    tt = dynInd(tt, idxTempFile, 2, yTempFile  );%1 since idx{1}==idx{2}
                    tt = resSub(tt, 2:3, N(2:3)) ;

                    for n=2:3
                       tt = ifftshiftGPU(tt,n);
                       tt=ifftGPU(tt,n)*N(n);
                       tt=fftshiftOperator(tt,2,0,n);
                    end
                    rec.y = cat(5, rec.y, tt);

                    %Reshape back into single k-space
                    tt = zerosL( resPop(ss{1}.pFlat, 2:3,[],2) );
                    tt = dynInd(tt,idxTempFile, 2, pTempFile  );%1 since idx{1}==idx{2}
                    tt = resSub(tt, 2:3, N(2:3));

                    for n=2:3
                       tt=ifftshiftGPU(tt,n);
                       tt=ifftGPU(tt,n)*N(n);
                       tt=fftshiftOperator(tt,2,0,n);
                    end
                    rec.PT.pSliceImage = cat(5, rec.PT.pSliceImage, tt);

                    %Compone PT time test
                    rec.PT.pTimeTest = cat(2, rec.PT.pTimeTest , pTimeTempFile);  
                end
                %Make filename
                if shift==0
                    suffExtra='';
                else
                    suffExtra=sprintf('_version%d',shift+1);
                end
                fileName = sprintf('Hybrid_%s_midres_%.0f_Int%d%s%s',resol,perc, useInt,suff,suffExtra);

                rec.Names.Name = fileName;

                %Save
                save(strcat(pathIn{1},fileName,'.mat'),'rec','-v7.3');
            end
        end
    end

end