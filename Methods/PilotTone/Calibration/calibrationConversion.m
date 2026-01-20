
function [AOut, TOut, TOutFit] = calibrationConversion( AIn, p, MTIn, NIn, MTOut, NOut, useOffset, deb)

%CALIBRATIONCONVERSION Converts the Pilot Tone calibration matrix from the logical coordinates (IJK) (units of voxels and radians around the centre), to
%   transformation parameters in another reference frame (defaults to the Right-Anterior-Superior (RAS) frame w.r.t. the scanner's iso-centre).
%   [ANEW, TNEW]=CALIBRATIONCONVERSION(AIN, P, MTIN, NIN, {MTOUT}, {NOUT}, {DI}, {BLOK})
%   * AIN is the current calibration matrix.
%   * P is the Pilot Tone (PT) signal.
%   * MTIN is the current s-form matrix.
%   * NIN is the current array size.
%   * {MTOUT} is the output s-form matrix for which to define the calibration matrix.
%   * {NOUT} is the output array size for which to define the calibration matrix.
%
%   Yannick Brackenier

if nargin < 5 || isempty(MTOut); MTOut = []; end
if nargin < 6 || isempty(NOut); NOut = []; end
if nargin < 7 || isempty(useOffset); useOffset = 0; end
if nargin < 8 || isempty(deb); deb = 0; end

%%% CREATE MOTION PARAMETERS FROM THE PILOT TONE SIGNAL
T = gather( pilotTonePrediction (AIn, p, 'PT2T','backward', size(p,1),size(AIn,1))  );


%%% CONVERT THE MOTION PARAMETERS
if isequal(MTIn,MTOut) && isequal(NIn,NOut) 
    fprintf('calibrationConversion:: Geometry the same so input parameters returned.\n')
    AOut=AIn;
    TOut=T;
    TOutFit=T;
else
    fprintf('calibrationConversion:: Converting calibration matrix from  MT_in\n')
    MTIn
    fprintf('To MT_out\n');
    MTOut
    TOut = transformationConversion_v2(T, MTIn, MTOut, NIn, NOut);
    
    %%% FIT CALIBRATION
    [AOut,~, TOutFit] = pilotToneCalibration(p, TOut, 'backward', useOffset);

    %%% PLOT
    if deb
        visMotionExt(TOut-TOutFit,[],[],0);
    end
end





