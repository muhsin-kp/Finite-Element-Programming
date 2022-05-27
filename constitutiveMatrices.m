function [AMat,BMat,DMat,SMat] = constitutiveMatrices(thickness,ks,pi,phi,orientationMat,laminateLayers,materialProperties)
    h=thickness;
    t=h/laminateLayers;

    %Assign Material Properties
  
        EX=materialProperties(1);EY=materialProperties(2);GXY=materialProperties(3);GYZ=materialProperties(4);GXZ=materialProperties(5);
        nyuXY=materialProperties(6);nyuYX=materialProperties(7);factor=1-nyuXY*nyuYX;
    
    % orienation of laminates

     if orientationMat==0
        orientationMat=[0,pi/2,0];
     end

    % Coordinates in thickness direction

    if rem(laminateLayers,2)==0
        for zz=1:laminateLayers+1
            z(zz)=(-laminateLayers*t/2)+(zz-1)*t;
        end
    elseif rem(laminateLayers,2)~=0
        for zz=1:laminateLayers+1
            z(zz)=(-laminateLayers*t/2)+(zz-1)*t;
        end
    end

    % Q Matrix

    q(1,1)=EX/factor;
    q(1,2)=nyuYX*EX/factor;
    q(2,1)=nyuXY*EY/factor;
    q(2,2)=EY/factor;
    q(3,3)=GXY;
    q(4,4)=ks*GYZ;
    q(5,5)=ks*GXZ;
    
    % Transformation matrix

    T=[cos(phi)^2,sin(phi)^2,-sin(2*phi),0,0;...
       sin(phi)^2,cos(phi)^2,sin(2*phi),0,0;...
       sin(phi)*cos(phi),-sin(phi)*cos(phi),cos(phi)^2-sin(phi)^2,0,0;...
       0,0,0,cos(phi),sin(phi);...
       0,0,0,-sin(phi),cos(phi)];
    
    % Q bar ,atrix

    qB=T*q*T.';
 
    % Substitute orientaion angle for phi

    for s=1:size(orientationMat,2)
        for i=1:5
            for j=1:5
               Q(i,j,s)=subs(qB(i,j,1),phi,orientationMat(s));
           end
       end
       Q=double(Q);
    end
    
% A,B,D,S matrix calculation

    Astiff(5,5)=0;Bstiff(5,5)=0;Dstiff(5,5)=0;Sstiff(5,5)=0;
    for k=1:size(orientationMat,2)
            for i=1:3
            for j=1:3
            Astiff(i,j)=Astiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
            Bstiff(i,j)=Bstiff(i,j)+Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
            Dstiff(i,j)=Dstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
            end
            end
    
                for i=4:5
                for j=4:5
                Sstiff(i,j)=Sstiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
                end
                end
    end

% A,B,D,S matrices

    AMat=Astiff(1:3,1:3);
    BMat=Bstiff(1:3,1:3);
    DMat=Dstiff(1:3,1:3);
    SMat=Sstiff(4:5,4:5);

end