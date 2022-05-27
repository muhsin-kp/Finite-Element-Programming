function [LengthL,LengthW,ThicknessZ,laminateLayers,ElementsX,ElementsY,Type,orientationMat,boundryCond,materialProperties,rho] = inputgui()
    %% 

    % Enter Geometry
    Geometry=str2double(inputdlg({'Plate Length','Plate Width','Total Plate Thickness','Number of laminas'},'Geometry Parametres'));
    
    %Enter Mesh Parameters
    Mesh=str2double(inputdlg({'Number of Divisions in X','Number of Divisions in X'},'Mesh Parametres'));
    
    %Enter Boundry Condition
    list1={'Clamped','Simply Supported'}
    [boundryCond,tf1]=listdlg('PromptString','Select Boundry Condition','ListString',list1,'SelectionMode','single');
   
    
    %Enter Type
    list2={'Cross-Ply','Hybrid'};
    [Type,tf2]=listdlg('PromptString','Select Type','ListString',list2,'SelectionMode','single');
 %% 
 
    
    %Enter Material Properties according to Type
    if tf2==1
        if Type==1   
            %orientationMat=inputdlg('Orientation Matrix (eg:[0,pi/2,0])','Orientation Matrix');
            mat=inputdlg({'EX','EY','GXY','GYZ','GXZ','NYUXY','NYUYX','Density'},'Enter material properties');
            materialProperties=str2double(mat);
            orientationMat=zeros(1,Geometry(4));
            for l=1:Geometry(4)
                text=strcat(num2str(l),' Layer Orientation Angle');
                orm=inputdlg(text,'Orientation Angle',[1 30],{'0'});
                orientationMat(1,l)=str2num(orm{1});
            end
        elseif Type==2
            materialProperties=zeros(Geometry(4),8);
            orientationMat=0;
            for l=1:Geometry(4)
                 mat=inputdlg({'EX','EY','GXY','GYZ','GXZ','NYUXY','NYUYX','Density'},sprintf('Enter %f Material Properties',num2str(l)));
                 materialProperties(l,:)=str2double(mat);
            end
        end
    end
    %% 

    LengthL=Geometry(1);
    LengthW=Geometry(2);
    ThicknessZ=Geometry(3);
    laminateLayers=Geometry(4);
    ElementsX=Mesh(1);
    ElementsY=Mesh(2);
    if boundryCond==1
        boundryCond='c-c-c-c';
    elseif boundryCond==2
        boundryCond='s-s-s-s';
    end
    if Type==1
        Type='C';
    elseif Type==2
        Type='H';
    end
    if isnan(materialProperties)
        materialProperties=[0,0,0,0,0,0,0,0];
        rho=1;
    else
        rho=materialProperties(8);
    end

end