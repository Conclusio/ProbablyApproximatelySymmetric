function [fullError,sampledError] = CalcFullAndSampledError3D(I1,I2,config,roiMask,fullOrSampled,epsilon)
%
% Calculates either full or sampled error or both depending on
% 'fullOrSampled' parameter. Parameter 'config' typically specifies only
% one config to test against. Full error means take the config and test it
% against all points (typically a lot) in the ROI. Sampled error means take
% the config and test against random sampled points (depending on epsilon).
% The error is calculated by first creating the affine transformation
% matrix from the config and then using this matrix against the points and
% measuring the distortion between the transformed points and the actual
% point at that location. The distortion is calculated per config.
%


computeFull = 0;
computeSampled = 0;

switch fullOrSampled
    case 'full'
        computeFull = 1;
        
    case 'sampled'
        computeSampled = 1;
        
    case 'both'
        computeFull = 1;
        computeSampled = 1;
end

% create matrix version of config
[h1,w1,d1] = size(I1); r1x = 0.5*(w1-1); r1y = 0.5*(h1-1); r1z = 0.5*(d1-1);
[h2,w2,d2] = size(I2); r2x = 0.5*(w2-1); r2y = 0.5*(h2-1); r2z = 0.5*(d2-1);
matrixConfig_mex = ...
    Configs2Affine_mex_3D(config',int32(h1),int32(w1),int32(d1),int32(h2),int32(w2),int32(d2),...
    int32(r1x),int32(r1y),int32(r1z),int32(r2x),int32(r2y),int32(r2z));

photometricInvariance = false;
roiIdxs = find(roiMask)';

% "full" roi sampling
if (computeFull)
    [ys,xs,zs] = ind2sub([h1,w1,d1],roiIdxs);
    
    fullError = EvaluateConfigs_mex_3D(I1,I2,matrixConfig_mex,int32(xs),int32(ys),int32(zs),int32(photometricInvariance));
else
    fullError = -1;
end

if (computeSampled)
    if (~exist('epsilon','var'))
        epsilon = 0.1; % a default value
    end
    numPoints = round(3/epsilon^2);
    idxs = randsample(roiIdxs,numPoints);
    [ys,xs,zs] = ind2sub([h1,w1,d1],idxs);
    sampledError = EvaluateConfigs_mex_3D(I1,I2,matrixConfig_mex,int32(xs),int32(ys),int32(zs),int32(photometricInvariance));
else
    sampledError = -1;
end