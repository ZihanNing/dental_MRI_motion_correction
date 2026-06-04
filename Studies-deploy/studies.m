% Auto-generated studies.m
% Case 1, file meas_MID03419_FID159770_PDwSPACE_DISORDER_ETL30_move

pathData{1}='./Studies-deploy/1/';
pathIn{1}=pathData{1};
pathRef{1}=strcat(pathData{1},'/Re-Se');
pathRemote{1}=pathData{1};

% ACQUISITIONS DATA
fileIn{1}{1}='meas_MID03419_FID159770_PDwSPACE_DISORDER_ETL30_move';

% REFERENCE DATA
refIn  = fillCell(fileIn, '');
refBIn = fillCell(fileIn, '');

% B0 and B1
B0In = fillCell(fileIn, '');
B1In = fillCell(fileIn, '');

% PILOT TONE FLAG
isPT     = fillCell(fileIn, 0);
facFOVTh = fillCell(fileIn, -.45);

% DATA CONVERSION SPECIFICATIONS
noiseFile      = fillCell(fileIn, refIn{1}{1});
supportReadout = fillCell(fileIn, []);
resRec         = fillCell(fileIn, []);
RDesired       = fillCell(fileIn, []);
