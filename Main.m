%FEM_COURSE PROJECT
%VIBRATION ANALYSIS OF LAMINATED PLATES
%DONE BY MUHSIN KP,VISHNU A
%% Pre-Processing

clc; clearvars; close all;ks=pi*pi/12;syms phi;

[LengthL,LengthW,thicknessZ,laminateLayers,ElementsX,ElementsY,laminateType,orientationMat,boundryCond,materialProperties,rho] = inputgui();


% Calculate constitutive matrix according type chosen by user, C=crossply
% H=Hybrid plates
      
    if laminateType=='C'
        [Amat,Bmat,Dmat,Smat]=constitutiveMatrices(thicknessZ,ks,pi,phi,orientationMat,laminateLayers,materialProperties);
    elseif laminateType=='H'
        [Amat,Bmat,Dmat,Smat] = constitutiveMatrices_sandwichMaterial(thicknessZ,ks,laminateLayers,materialProperties);
        % average density of laminates
        rho=sum(materialProperties(:,8),1)/laminateLayers;
    end


%Generating Mesh

tic
[totalElements,totalNodes,nodalCoordinateMat,elementConnectivityMat] = generateMesh(LengthL,LengthW,...
                                         ElementsX,ElementsY);

GlobalDof=5*totalNodes; 

%Generating Stiffness Matrix

[stiffness]=calculateStiffnessMatrix(GlobalDof,totalElements,...
    elementConnectivityMat,totalNodes,nodalCoordinateMat,Amat,Bmat,Dmat,Smat);

%Generating Mass Matrix


[mass]=calculateMassMatrix(GlobalDof,totalElements,...
    elementConnectivityMat,totalNodes,nodalCoordinateMat,rho,thicknessZ);


% default BC

if isempty(boundryCond)
    boundryCond='c-c-c-c';
end

% enforcing boundary conditions
[unConstrainedDOFs,constrainedNodesinW]=...
    boundryConditions(boundryCond,GlobalDof,nodalCoordinateMat,totalNodes);

%% Processing


[eigVectors,eigValues] = eig(stiffness(unConstrainedDOFs,unConstrainedDOFs),...
    mass(unConstrainedDOFs,unConstrainedDOFs)); 

%% Post-Processing

  reqModes=6;  
  plotModeShape(eigValues,eigVectors,thicknessZ,reqModes,totalNodes,LengthL,LengthW,ElementsX,ElementsY,constrainedNodesinW,materialProperties)
time=toc;
