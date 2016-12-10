function A = CreateAffineTransformation3D(trans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tx =  trans(1);
% ty =  trans(2);
% r2 =  trans(3);
% sx =  trans(4);
% sy =  trans(5);
% r1 =  trans(6);
%
%
% % syms tx ty r2 sx sy r1
% % syms R1 T2 S A
%
% R1 = [ cos(r1) , -sin(r1) ; sin(r1) , cos(r1) ];
% R2 = [ cos(r2) , -sin(r2) ; sin(r2) , cos(r2) ];
% S =  [ sx      , 0        ; 0       , sy      ];
% A =  [0 0 tx ; 0 0 ty ; 0 0 1];
% A(1:2,1:2) = R1*S*R2;


% A =     [ sx*cos(r1)*cos(r2) - sy*sin(r1)*sin(r2), - sx*cos(r1)*sin(r2) - sy*cos(r2)*sin(r1), tx]
%         [ sx*cos(r2)*sin(r1) + sy*cos(r1)*sin(r2),   sy*cos(r1)*cos(r2) - sx*sin(r1)*sin(r2), ty]
%         [                                       0,                                         0,  1]

% disp('done');

tx =  trans(1);
ty =  trans(2);
tz =  trans(3);
s =  trans(4);
latitude =  trans(5);
longitude =  trans(6);
roll = trans(7);


% these are the coordinates of the unit sphere point, around which a rotation by 'roll' is performed
ux = sin(latitude)*cos(longitude);
uy = sin(latitude)*sin(longitude);
uz = cos(latitude);

sinRoll = sin(roll);
cosRoll = cos(roll);

% and the matrix
A = zeros(4);

if (s < 0)    
    % the case of reflection w.r.t. the normal [ux,uy,uz]:
    % A =   [1-2uxux  -2uxuy  -2uxuz ;
    %       [-2uxuy   1-2uyuy -2uyuz ;
    %       [-2uxuz   -2uyuz  1-2uzuz]
    A(1,1) = 1 - 2*ux*ux;
    A(1,2) = -2*ux*uy;
    A(1,3) = -2*ux*uz;
    A(1,4) = tx;
    A(2,1) = -2*uy*ux;
    A(2,2) = 1-2*uy*uy;
    A(2,3) = -2*uy*uz;
    A(2,4) = ty;
    A(3,1) = -2*uz*ux;
    A(3,2) = -2*uz*uy;
    A(3,3) = 1-2*uz*uz;
    A(3,4) = tz;
else    
    A(1,1) = s*(cosRoll+ux*ux*(1-cosRoll));
    A(1,2) = s*(ux*uy*(1-cosRoll)-uz*sinRoll);
    A(1,3) = s*(ux*uz*(1-cosRoll)+uy*sinRoll);
    A(1,4) = tx;
    A(2,1) = s*(uy*ux*(1-cosRoll)+uz*sinRoll);
    A(2,2) = s*(cosRoll+uy*uy*(1-cosRoll));
    A(2,3) = s*(uy*uz*(1-cosRoll)-ux*sinRoll);
    A(2,4) = ty;
    A(3,1) = s*(uz*ux*(1-cosRoll)-uy*sinRoll);
    A(3,2) = s*(uz*uy*(1-cosRoll)+ux*sinRoll);
    A(3,3) = s*(cosRoll+uz*uz*(1-cosRoll));
    A(3,4) = tz;
end

A(4,1) = 0;
A(4,2) = 0;
A(4,3) = 0;
A(4,4) = 1;








