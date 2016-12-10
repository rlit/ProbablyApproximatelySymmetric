function thresh = GetThreshPerDelta3D(delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on experimentally drawn % p = 0.1341    0.0278

p = [0.2, 0.0278];
safety = p(2);%0.03;
thresh = p(1)*delta+p(2)-safety;

