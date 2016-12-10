function dist1 = VisualizeSamplingGrayscale3D(resultName,shape1,shape2,...
    numPoints,a,tx,ty,tz,xs,ys,zs,config,optA,...
    sampledError,fullError,rotationFold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('optA','var'))
    optA = [];
end

if isstruct(shape1)
    plotMesh = true;
    shape1 = BringMeshToPixelCoords(shape1);
    shape2 = BringMeshToPixelCoords(shape2);
    I1 = shape1.G;
    I2 = shape2.G;
    
    if isfield(shape1,'plotBoth') && isequal(shape1.plotBoth,1)
        plotBothShapes = true;
    else
        plotBothShapes = false;
    end
    
    
    % shape 2 boulding box
    bb2.xMin = min(shape2.X);
    bb2.xMax = max(shape2.X);
    bb2.yMin = min(shape2.Y);
    bb2.yMax = max(shape2.Y);
    bb2.zMin = min(shape2.Z);
    bb2.zMax = max(shape2.Z);
else
    plotMesh = false;
    I1 = shape1;
    I2 = shape2;
    
    % shape 2 boulding box
    szMax = max(size(I1),size(I2));
    bb2.xMin = 1;
    bb2.yMin = 1;
    bb2.zMin = 1;
    bb2.xMax = szMax(2);
    bb2.yMax = szMax(1);
    bb2.zMax = szMax(3);
end


a3x3 = a(1:3,1:3);

%% error  calculation
if ~(exist('sampledError','var') && exist('fullError','var'))
    [dist,dist1] = CalculateError(I1,I2,numPoints,a3x3,tx,ty,tz,xs,ys,zs);
end

%% the figure

figure; hold on; axis equal;
facealpha = 0.5;

%color_I1 = 'r';
color_I2 = 'g';
color_I1_trans = 'b';
color_OPT = 'y';



if plotMesh

    [shape1_trans] = TransformMesh(shape1,a3x3,tx,ty,tz);
    if plotBothShapes
        plotFun = @(s,c,a,n)trisurf(s.TRIV, s.Y, s.X, s.Z,'facecolor',c,'edgecolor','none','facealpha',a,'displayname',n);
        %plotFun(shape1,color_I1,facealpha)
        plotFun(shape2,      color_I2,      facealpha,'shape_2')
        plotFun(shape1_trans,color_I1_trans,facealpha,'shape_1_trans')
    else
        dtVals  = GetMeshErrorVals(shape1_trans,shape2);
        pathcH = trisurf(shape2.TRIV, shape2.Y, shape2.X, shape2.Z,dtVals,'displayname','shape_error');
        set(pathcH,'edgecolor','none')
        shape2_radius = (size(shape2.G,1)-1)/2; % <-- assumes boxified shape
        set(gca,'cLim',[0 shape2_radius/10]);
        shading interp
    end

    if (~isempty(optA))
        error('not implemented for mesh')
    end

else
    % perform a backwards mapping (NN) + create isosurface

    [targetMatBack,targetMatBack_OPT] = DoBackwardMap(I1,I2,a3x3,tx,ty,tz,optA);
    plotFun = @(I,c,a)patch(isosurface(I,0.5),'facecolor',c,'edgecolor','none','facealpha',a);

    %plotFun(I1,color_I1,facealpha);
    plotFun(I2,color_I2,facealpha);
    plotFun(targetMatBack,color_I1_trans,facealpha);

    if (~isempty(optA))
    plotFun(targetMatBack_OPT,color_OPT,facealpha);
    end

end

%% plot the config
titlePrefix = DepictConfig(config,rotationFold,size(I1));

%%
view(3);
camlight
% % set the view from the rotation axis
% view([axis_x,axis_y,axis_z])

if ~(exist('sampledError','var') && exist('fullError','var'))
    title([resultName ' ' titlePrefix '  with sampled dist = ' num2str(dist)])
else
    title([resultName ' ' titlePrefix '  with sampled dist = ' num2str(sampledError) ' , full dist = ' num2str(fullError)])
end

axis([bb2.yMin bb2.yMax bb2.xMin bb2.xMax bb2.zMin bb2.zMax]);

xlabel X; ylabel Y; zlabel Z;



function [dist,dist1] = CalculateError(I1,I2,numPoints,a3x3,tx,ty,tz,xs,ys,zs)
%% sizes
[h1,w1,d1] = size(I1);
[h2,w2,d2] = size(I2);

r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r1z = 0.5*(d1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);
r2z = 0.5*(d2-1);

%% sampled source and target points
targetPoints = a3x3*[xs-(r1x+1);ys-(r1y+1);zs-(r1z+1)];

txs = round(targetPoints(1,:) + (r2x+1)  + tx);
tys = round(targetPoints(2,:) + (r2y+1)  + ty);
tzs = round(targetPoints(3,:) + (r2z+1)  + tz);
insideInds = find(tzs>0 & tzs<(d2+1) & txs>0 & txs<(w2+1) & tys>0 & tys<(h2+1));

%% ones of original shape
one_inds = find(I1);
[one_ys,one_xs,one_zs] = ind2sub(size(I1),one_inds');

ones_targetPoints = a3x3*[one_xs-(r1x+1);one_ys-(r1y+1);one_zs-(r1z+1)];

ones_txs = BoundBy(round(ones_targetPoints(1,:) + (r2x+1)  + tx),1,w2);
ones_tys = BoundBy(round(ones_targetPoints(2,:) + (r2y+1)  + ty),1,h2);
ones_tzs = BoundBy(round(ones_targetPoints(3,:) + (r2z+1)  + tz),1,d2);

targetMat = zeros(size(I2));
target_one_inds = sub2ind(size(I2),ones_tys,ones_txs,ones_tzs);
targetMat(target_one_inds) = 1;


%% scores
sourceInds = sub2ind([h1,w1,d1],ys(insideInds),xs(insideInds),zs(insideInds));
targetInds = sub2ind([h2,w2,d2],tys(insideInds),txs(insideInds),tzs(insideInds));
% badfitinds = insideInds(I1(sourceInds) ~= I2(targetInds));
% outsideinds = setdiff(1:numPoints,insideInds);
% badMatches = length(badfitinds);
badMatches =  sum(abs(I1(sourceInds) - I2(targetInds)));
outOfBounds = numPoints - length(insideInds);
dist = (badMatches + outOfBounds)/numPoints;
dist1 = badMatches/numPoints;

function titlePrefix = DepictConfig(config,rotationFold,imSize)
r1x = 0.5*(imSize(2)-1);
r1y = 0.5*(imSize(1)-1);
r1z = 0.5*(imSize(3)-1);

s   =  config(4);
lat = -config(5)+pi/2;
lon =  config(6);

[axis_x,axis_y,axis_z] = sph2cart(lon,lat,r1z*1.4);
if (s ~=-1)
    % rotation
    if (exist('rotationFold','var') && rotationFold > 1)
        axisColor = 'r';
        titlePrefix = sprintf('%d-FOLD rotation ',rotationFold);
    elseif (exist('rotationFold','var') && rotationFold < 2) % revolutional
        axisColor = 'm';
        titlePrefix = sprintf('REVOLUTION rotation ');
    else
        axisColor = 'r';
        titlePrefix = 'REGULAR rotation';
    end
    
    axisPoint1 = [r1x r1y r1z] - [axis_x,axis_y,axis_z];
    axisPoint2 = [r1x r1y r1z] + [axis_x,axis_y,axis_z];
    R_ = r1x / 70;
    N_ = 20;
    [CylinderHandle] = Cylinder(axisPoint1,axisPoint2,R_,N_,axisColor);
    set(CylinderHandle,'displayname','rotation_axis')
else
    % reflection

    titlePrefix = 'regular reflection';

    % 1st basis vector for the reflection plane: select smallest coordinate
    % 2nd basis vector for the reflection plane: orthogonal to latter & the reflection axis
    [~,idx] = min(abs([axis_x,axis_y,axis_z]));
    tempVec = zeros(3,1);tempVec(idx) = 1;
    ortho1 = cross([axis_x,axis_y,axis_z],tempVec); ortho1 = ortho1/norm(ortho1)*r1z*1;
    ortho2 = cross([axis_x,axis_y,axis_z],ortho1);  ortho2 = ortho2/norm(ortho2)*r1z*1;
    
    [axispoint1,axisPoint2] = meshgrid(-1:.1:1);
    GRD = [axispoint1(:),axisPoint2(:)];
    
    fv.faces    = delaunay(axispoint1(:),axisPoint2(:));
    fv.vertices = [...
        r1x+1 + ortho1(1)*GRD(:,1) + ortho2(1)*GRD(:,2),...
        r1y+1 + ortho1(2)*GRD(:,1) + ortho2(2)*GRD(:,2),...
        r1z+1 + ortho1(3)*GRD(:,1) + ortho2(3)*GRD(:,2)];
    p = patch(fv);
    facealpha = 0.5;
    set(p,'faceColor','r','facealpha',facealpha,'edgeColor','none','displayname','reflection_plane')

end

function [targetMatBack,targetMatBack_OPT] = DoBackwardMap(I1,I2,a3x3,tx,ty,tz,optA)
%% sizes
[h1,w1,d1] = size(I1);
[h2,w2,d2] = size(I2);

r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r1z = 0.5*(d1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);
r2z = 0.5*(d2-1);

%% perform a backwards mapping (NN)
[Xs,Ys,Zs] = meshgrid(1:w2,1:h2,1:d2);
Xs = Xs - (r2x+1) - tx;
Ys = Ys - (r2y+1) - ty;
Zs = Zs - (r2z+1) - tz;
coords = [Xs(:)';Ys(:)';Zs(:)'];
sourceVoxs = a3x3\coords;
clear coords Xs Ys Zs
sv_X = sourceVoxs(1,:)+(r1x+1);
sv_Y = sourceVoxs(2,:)+(r1y+1);
sv_Z = sourceVoxs(3,:)+(r1z+1);
reshaped_sv_X = reshape(sv_X,size(I2));
reshaped_sv_Y = reshape(sv_Y,size(I2));
reshaped_sv_Z = reshape(sv_Z,size(I2));

clear sourceVoxs sv_X sv_Y sv_Z

targetMatBack = interp3(I1,reshaped_sv_X,reshaped_sv_Y,reshaped_sv_Z,'nearest',0);

% and using the OPT:
if (~isempty(optA))
    [Xs,Ys,Zs] = meshgrid(1:w2,1:h2,1:d2);
    a3x3_OPT = optA(1:3,1:3);
    tx_OPT   = optA(1,4);
    ty_OPT   = optA(2,4);
    tz_OPT   = optA(3,4);
    Xs = Xs - (r2x+1) - tx_OPT;
    Ys = Ys - (r2y+1) - ty_OPT;
    Zs = Zs - (r2z+1) - tz_OPT;
    coords = [Xs(:)';Ys(:)';Zs(:)'];
    sourceVoxs = a3x3_OPT\coords;
    clear coords Xs Ys Zs
    sv_X = sourceVoxs(1,:)+(r1x+1);
    sv_Y = sourceVoxs(2,:)+(r1y+1);
    sv_Z = sourceVoxs(3,:)+(r1z+1);
    reshaped_sv_X = reshape(sv_X,size(I2));
    reshaped_sv_Y = reshape(sv_Y,size(I2));
    reshaped_sv_Z = reshape(sv_Z,size(I2));

    clear sourceVoxs sv_X sv_Y sv_Z

    targetMatBack_OPT = interp3(I1,reshaped_sv_X,reshaped_sv_Y,reshaped_sv_Z);
else
    targetMatBack_OPT = [];
end

function [shape1_trans] = TransformMesh(shape1,a3x3,tx,ty,tz)
s = shape1;

[h1,w1,d1] = size(s.G);
r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r1z = 0.5*(d1-1);

s.X = s.X - (r1y+1) - ty;
s.Y = s.Y - (r1x+1) - tx;
s.Z = s.Z - (r1z+1) - tz;
coords = [s.Y s.X s.Z ]';
sourceVoxs = a3x3\coords;
s.X = sourceVoxs(2,:)+(r1y+1);
s.Y = sourceVoxs(1,:)+(r1x+1);
s.Z = sourceVoxs(3,:)+(r1z+1);


[s.G] = DoBackwardMap(s.G,s.G,a3x3,tx,ty,tz,[]);
shape1_trans = s;


% figure(2345);clf
% hold on
% % plotFun = @(s,c,a)trisurf(s.TRIV, s.Y, s.X, s.Z,'facecolor',c,'edgecolor','none','facealpha',a);
% % plotFun(shape1,'r',.5)
% % plotFun(s,'g',.5)
% plotFun = @(I,c,a)patch(isosurface(I,0.5),'facecolor',c,'edgecolor','none','facealpha',a);
% plotFun(shape1.G,'r',.5)
% plotFun(s.G,'g',.5)
%
% view(3)
% axis equal on
% grid on




function shapeOut = BringMeshToPixelCoords(shapeIn)
shapeOut = shapeIn;
Img1 = shapeOut.G;
%% clipping
idx = find(Img1);
[Ys,Xs,Zs] = ind2sub(size(Img1),idx);
minShift = [min(Ys), min(Xs), min(Zs)];
I_clip = Img1(min(Ys):max(Ys) , min(Xs):max(Xs) , min(Zs):max(Zs));

%% COG % radius
idx = find(I_clip);
[Ys,Xs,Zs] = ind2sub(size(I_clip),idx);
COG = round(mean([Ys,Xs,Zs]));
diffs = bsxfun(@minus,[Ys,Xs,Zs],COG)';
dists = sqrt(sum(diffs.^2));
radius = ceil(max(dists));

sz = size(I_clip);
rngFun = @(i)(radius+1-(COG(i)-1)) : (radius+1-(COG(i)-1) + sz(i)-1);
yRng = rngFun(1);
xRng = rngFun(2);
zRng = rngFun(3);
shapeOut.G = zeros([2*radius+1,2*radius+1,2*radius+1]);
shapeOut.G(yRng,xRng,zRng) = I_clip;

%%

shapeOut.X = shapeOut.X / shapeOut.d;
shapeOut.Y = shapeOut.Y / shapeOut.d;
shapeOut.Z = shapeOut.Z / shapeOut.d;

shapeOut.X = shapeOut.X - min(shapeOut.X) + 1 + yRng(1) - minShift(1);
shapeOut.Y = shapeOut.Y - min(shapeOut.Y) + 1 + xRng(1) - minShift(2);
shapeOut.Z = shapeOut.Z - min(shapeOut.Z) + 1 + zRng(1) - minShift(3);
%%
% figure(2345);clf
% hold on
% plotFun = @(s,c,a)trisurf(s.TRIV, s.Y, s.X, s.Z,'facecolor',c,'edgecolor','none','facealpha',a);
% plotFun(shapeIn,'r',.5)
% plotFun = @(I,c,a)patch(isosurface(I,0.5),'facecolor',c,'edgecolor','none','facealpha',a);
% plotFun(shapeOut.G,'g',.5);
% view(3)
% axis equal

function errVals = GetMeshErrorVals(shape1,shape2)
shape2_DT = double(bwdist(bwperim(shape1.G)));
errVals   = interp3(shape2_DT,shape2.Y,shape2.X,shape2.Z);

%replace NaN with maximal value
errVals(isnan(errVals)) = max(errVals(~isnan(errVals)));

%truncate at 99%
pct99 = prctile(errVals,99);
errVals(errVals>pct99) = pct99;


