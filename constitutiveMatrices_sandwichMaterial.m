function [AMat,BMat,DMat,SMat] = constitutiveMatrices_sandwichMaterial(thicknessZ,ks,laminateLayers,materialProperties)

   h=thicknessZ;
   t=h/laminateLayers;


    if rem(laminateLayers,2)==0
        for zz=1:laminateLayers+1
            z(zz)=(-laminateLayers*t/2)+(zz-1)*t;
        end
    elseif rem(laminateLayers,2)~=0
        for zz=1:laminateLayers+1
            z(zz)=(-laminateLayers*t/2)+(zz-1)*t;
        end
    end

    for l=1:laminateLayers
        factor=1-materialProperties(l,7)*materialProperties(l,6);
        q(1,1,l)=materialProperties(l,1)/factor;
        q(1,2,l)=materialProperties(l,7)*materialProperties(l,1)/factor;
        q(2,1,l)=materialProperties(l,6)*materialProperties(l,2)/factor;
        q(2,2,l)=materialProperties(l,2)/factor;
        q(3,3,l)=materialProperties(l,3);
        q(4,4,l)=ks*materialProperties(l,4);
        q(5,5,l)=ks*materialProperties(l,5);
    end    
    Q=q;
    
    %A,B,D,S Matrix calculation

    Astiff(5,5)=0;Bstiff(5,5)=0;Dstiff(5,5)=0;Sstiff(5,5)=0;
    
            for i=1:3
                for j=1:3
                    for k=1:laminateLayers
                        Astiff(i,j)=Astiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
                        Bstiff(i,j)=Bstiff(i,j)+Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
                        Dstiff(i,j)=Dstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
                    end
               end
            end
    
                for i=4:5
                    for j=4:5
                        for k=1:laminateLayers
                            Sstiff(i,j)=Sstiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
                        end
                    end
                end
    

    AMat=Astiff(1:3,1:3);
    BMat=Bstiff(1:3,1:3);
    DMat=Dstiff(1:3,1:3);
    SMat=Sstiff(4:5,4:5);
