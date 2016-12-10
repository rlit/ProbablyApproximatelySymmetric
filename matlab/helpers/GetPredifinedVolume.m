function [I1,I1_radius] = GetPredifinedVolume(shape)
maxDim = 151;
I1_temp = zeros(maxDim * [1 1 1]);

icosahedron = 0;
dodecahedron = 0;

if (exist('shape','var'))
    switch shape
        case 'icosahedron'
            icosahedron = 1;
        case 'dodecahedron'
            dodecahedron = 1;
    end
end
if (icosahedron||dodecahedron)

    %%
    phi = (1 + sqrt(5)) / 2; % the golden ratio

    if (icosahedron)
        xs = [0     0    0    0   phi phi -phi -phi    1    -1    1   -1];
        ys = [1    -1    1   -1   0     0    0    0  phi phi -phi -phi];
        zs = [phi phi -phi -phi   1    -1    1   -1    0   0    0    0];
        vertexDist = 2;
    else
        xs = [1     1    1    1  -1  -1   -1   -1];
        ys = [1     1   -1   -1   1   1   -1   -1];
        zs = [1    -1    1   -1   1  -1    1   -1];

        xs = [xs [    0       0     0      0     phi     phi   -phi    -phi  1/phi  -1/phi  1/phi -1/phi]];
        ys = [ys [1/phi  -1/phi 1/phi -1/phi       0       0      0       0    phi     phi   -phi   -phi]];
        zs = [zs [phi       phi  -phi   -phi   1/phi  -1/phi  1/phi  -1/phi      0       0      0      0]];
        vertexDist = 2/phi;
    end

    fact = 30;
    xs = fact*(2+xs);
    ys = fact*(2+ys);
    zs = fact*(2+zs);
    shape = zeros(fact*[4,4,4]);
    vertices = [xs' ys' zs'];


    for source = 1 : length(xs)
        for target = source + 1 : length(xs)
            sourceVec = vertices(source,:);
            targetVec = vertices(target,:);
            dist = sqrt(sum2((sourceVec-targetVec).^2));
            if dist < fact*(1.01*vertexDist)
                [X,Y,Z] = bresenham_line3d(sourceVec, targetVec);
                shape(sub2ind(size(shape),round(Y),round(X),round(Z))) = 1;
            end
        end
    end

    %%
    rad = 4;
    [xxs,yys,zzs] = meshgrid(-rad:rad,-rad:rad,-rad:rad);
    dists = sqrt(xxs.^2+yys.^2+zzs.^2);
    SE = dists <= rad;

    dshape = imdilate(shape,SE);
    % dshape = ordfilt3D(dshape,14);

    % display
    fv1 = isosurface(dshape,0.5);
    figure; hold on; axis equal;
    patch(fv1,'facecolor','r','edgecolor','none');
    view(3); camlight

    I1_temp = dshape;
    %%
else

    %

    % example 3: octagram
    cntr   = ceil(maxDim/2);
    xRange = (-50:50) + cntr;
    yRange = (-50:50) + cntr;
    zRange = (-15:15) + cntr;
    I1_temp(yRange,xRange,zRange) = 1;

    % %add leg
    %legWidth = 9;
    %tmp = (1:legWidth)-1;
    %I1_temp(min(yRange)+tmp,min(xRange)+tmp,yRange)=1;

    %I1_temp = imrotate(I1_temp,90,'crop');
    I2_temp = imrotate(I1_temp,45,'crop');
    I1_temp = I1_temp | I2_temp;

end


[I1] = BoxifyShape(I1_temp);


% rotate by some angle:
latitude = pi*0.85;
longitude = pi/1.53;
roll = pi/7;

I1 = RotateByAngle(I1,latitude,longitude,roll);

% re-boxify
[I1,I1_radius] = BoxifyShape(I1);

maxDim = max(size(I1));
% show input
fv1 = isosurface(I1,0.5);
figure; hold on; axis equal;
patch(fv1,'facecolor','r','edgecolor','none');
view(3); camlight

axis([1,size(I1,1),1,size(I1,2),1,size(I1,3)]);