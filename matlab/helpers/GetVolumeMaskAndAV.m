function [Img,roiMask,bounds,steps,ATV,times] = GetVolumeMaskAndAV(K,I1,params,I1_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns
% - Img ... the distance field or the given image I1 (if K <= 0)
% - roiMask ... the region of interest according to roiMaskType of the parameters
% - bounds ... yet unknown
% - steps ... yet unknown
% - ATV ... some kind of average
% - times ... debug info on how long specific code parts ran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


workOnDTs = K>0;
tStart = tic;
if workOnDTs
    % Note: DT .. Distance transform
    % create full DT volume
    I1dt_full = double(bwdist(bwperim(I1))); % create distance field
    I1dt_full(I1==0) = -I1dt_full(I1==0); % outer parts have negative distance (?)
    I1dt_full = round(I1dt_full);

    % truncate at k
    I1dt = I1dt_full;
    clear I1dt_full
    I1dt(abs(I1dt)>K) = K * sign(I1dt(abs(I1dt)>K)); % truncated to +/- K

    % normalize values to [0,1]
    I1dt = I1dt - min(I1dt(:));
    I1dt = I1dt / max(I1dt(:));
    Img = I1dt;
    clear I1dt
else
    Img = I1;
end
times.tsdf = toc(tStart); % measure run-time of truncated signed distance field creation
%% roi mask
switch params.roiMaskType
    case 'sphere'
        idx = 1:numel(I1);
        [Ys,Xs,Zs] = ind2sub(size(I1),idx');
        maxDim = size(I1,1);
        rad = (maxDim-1)/2;
        cntr = [rad rad rad]+1;
        diffs = bsxfun(@minus,[Ys,Xs,Zs],cntr)';
        dists = sqrt(sum(diffs.^2)); % euclidean distances from center point, which is the center of the whole grid
        dists = reshape(dists,size(I1));
        roiMask = dists < rad; % all points which are in a sphere around the center point with radius of half the length of the maximum dimension
    case 'interior'
        roiMask = Img>0;
    case 'interiorAndSphereIntersect'
        idx = 1:numel(I1);
        [Ys,Xs,Zs] = ind2sub(size(I1),idx');
        clear idx
        maxDim = size(I1,1);
        rad = (maxDim-1)/2;
        cntr = [rad rad rad]+1;
        diffs = bsxfun(@minus,[Ys,Xs,Zs],cntr)';
        dists = sqrt(sum(diffs.^2));
        dists = reshape(dists,size(I1));
        roiMask = dists < rad;
        roiMask = roiMask & (Img>0); % the intersection
    otherwise
        error('bad roiMaskType');
end

%% determine search limits and initial grid
tStart = tic;
[bounds,steps,ATV] = DetermineGridBoundsAndSteps(Img,I1_radius,params.delta,roiMask,params.translationOffset);
times.tv = toc(tStart);

return
