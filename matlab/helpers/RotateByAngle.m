function I1_rot = RotateByAngle(I1,latitude,longitude,roll)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I2 = I1;
[w2,h2,d2] = size(I2);

I1_rot = zeros(size(I2));
[w1,h1,d1] = size(I1_rot);


rad1X = (w1+1)/2;
rad1Y = (h1+1)/2;
rad1Z = (d1+1)/2;
rad2X = (w2+1)/2;
rad2Y = (h2+1)/2;
rad2Z = (d2+1)/2;

% params
tx = 0; ty = 0; tz = 0;
scale = 1;
if (~exist('latitude','var'))
    latitude = pi/4;
    longitude = pi/3;
    roll = pi/7;
end

A = CreateAffineTransformation3D([tx,ty,tz,scale,latitude,longitude,roll]);
a3x3 = A(1:3,1:3);

[Xs,Ys,Zs] = meshgrid(1:w1,1:h1,1:d1);
Xs = Xs - (rad1X+1);
Ys = Ys - (rad1Y+1);
Zs = Zs - (rad1Z+1);
coords = [Xs(:)';Ys(:)';Zs(:)'];
targetVoxelss = a3x3*coords;
tv_X = targetVoxelss(1,:)+tx+(rad2X+1);
tv_Y = targetVoxelss(2,:)+ty+(rad2Y+1);
tv_Z = targetVoxelss(3,:)+tz+(rad2Z+1);
reshaped_tv_X = reshape(tv_X,[h1,w1,d1]);
reshaped_tv_Y = reshape(tv_Y,[h1,w1,d1]);
reshaped_tv_Z = reshape(tv_Z,[h1,w1,d1]);

I1_rot = interp3(I2,reshaped_tv_X,reshaped_tv_Y,reshaped_tv_Z,'nearest',0);
