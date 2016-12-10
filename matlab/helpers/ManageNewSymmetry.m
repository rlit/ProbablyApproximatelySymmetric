function detectedSymmetries = ManageNewSymmetry(detectedSymmetries,newSymmetry,Img,roiMask,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(detectedSymmetries)
    detectedSymmetries = newSymmetry;
else
    detectedSymmetries(end+1) = newSymmetry;
end

% referenceScore = newSymmetry.fullError;

if (newSymmetry.isReflection)
    if (params.expandDetectedSymmetry)
        
        % check if there was a previous reflection - for revolutionary detection
        for i = 1 : length(detectedSymmetries)-1
            prevSymmetry = detectedSymmetries(i);
            if (~prevSymmetry.isReflection)
                continue
            end
            
            % create suspected rotation axis
            axis1 = Config2Cart(newSymmetry.config);
            axis2 = Config2Cart(prevSymmetry.config);
            newAxis = cross(axis1,axis2);
            [lat,lon] = LatLongFromCart(newAxis);
            
            % check if this 'potential' rotation axis is already in list
            tmpConfig = [0 0 0 1 lat lon -pi];
            farIdxs = FindConfigsFarFromOnesInList(tmpConfig,detectedSymmetries,params.minAngleBetweenSymmetries);
            if (isempty(farIdxs))
                continue;
            end
            
            % run rotation fold/revolution test
            
            [rotationFold,representativeConfig,representativeScore] = DetectRotationAxisFold(lat,lon,Img,roiMask,params);
            newRotationaSym = detectedSymmetries(end);
            newRotationaSym.config = representativeConfig;
            newRotationaSym.isReflection = 0;
            newRotationaSym.rotationFold = rotationFold;
            newRotationaSym.fullError    = representativeScore;
            newRotationaSym.sampledError = representativeScore;
            newRotationaSym.binFullError = NaN;
            newRotationaSym.binSampled   = NaN;
            
            if rotationFold > 1
                % add rotational symmetry
                detectedSymmetries(end+1) = newRotationaSym; %#ok<AGROW>
                
            else
                % delete the first reflection
                detectedSymmetries(i) = [];
                
                % replace 2nd reflection with "new" revolution
                detectedSymmetries(end) = newRotationaSym;
            end
            
        end
    end
    
else
    lat = newSymmetry.config(5);
    lon = newSymmetry.config(6);
    [rotationFold,representativeConfig,representativeScore] = ...
        DetectRotationAxisFold(lat,lon,Img,roiMask,params);
    
    detectedSymmetries(end).rotationFold = rotationFold;
    detectedSymmetries(end).fullError    = representativeScore;
    detectedSymmetries(end).sampledError = representativeScore;
    detectedSymmetries(end).config = representativeConfig;
    
end



function res_axis = Config2Cart(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = -config(5)+pi/2;
lon = config(6);

[axis_x,axis_y,axis_z] = sph2cart(lon,lat,1);
res_axis = [axis_x,axis_y,axis_z];


function [lat,long] = LatLongFromCart(axis_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[long,lat] = cart2sph(axis_in(1),axis_in(2),axis_in(3));
lat = -lat+pi/2;






