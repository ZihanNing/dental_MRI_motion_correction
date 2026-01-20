
clc; close all;
clear all
cd /home/ybr19/Software/Registration/elastixReg/
addpath(genpath('.'))
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))

[x1 , MS1, MT1 ] = readNII('x1',{''});
[x2 , MS2, MT2 ] = readNII('x2',{''});

x1 = x1{1};x2 = x2{1};
MS1 = MS1{1};MS2 = MS2{1};
MT1 = MT1{1};MT2 = MT2{1};

params=[];
params.HistBins='32';
params.NumRes='3';
params.NumIter='200';
params.SpatSam='1000';
params.regType = 'affine';
params.GridSpac='5';%Only used for bspline registration

[xReg, paramReg] = elastixRegistration( x1, x2, MT1, MT2, params);

y=[];y{1}=cat(4, x1, x2, xReg);
MS=[];MS{1}=MS1;
MT=[];MT{1}=MT1;
writeNII('x',{'Reg'}, y, MS, MT); 
