function [bounds,steps] = GenerateGrid3D(delta,minTx,maxTx,minTy,maxTy,minTz,maxTz,minR,maxR,minS,maxS,...
    roiAverageTV,roiRadius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% search bounds
bounds_tx = [minTx,maxTx];
bounds_ty = [minTy,maxTy];
bounds_tz = [minTz,maxTz];
bounds_r = [minR,maxR];
bounds_s = [minS,maxS]; % [1/c,c]%

bounds.tx = bounds_tx;
bounds.ty = bounds_ty;
bounds.tz = bounds_tz;
bounds.r = bounds_r;
bounds.s = bounds_s;

%% step sizes

MaximalVoxelMotion = delta/roiAverageTV;
step_tx = MaximalVoxelMotion/sqrt(3);
step_ty = MaximalVoxelMotion/sqrt(3);
step_tz = MaximalVoxelMotion/sqrt(3);

step_r = MaximalVoxelMotion/roiRadius;

step_s = MaximalVoxelMotion/roiRadius; % TODO - verify this

steps.tx = step_tx;
steps.ty = step_ty;
steps.tz = step_tz;
steps.r = step_r;
steps.s = step_s;
