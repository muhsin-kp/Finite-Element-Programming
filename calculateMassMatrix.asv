function [M]=calculateMassMatrix(GlobalDof,totalElements,elementConnectivityMat,totalNodes,nodalCoordinateMat,rho,thickness)


%Moment of inertia
I=thickness^3/12;

%Initiate Global Stiffness Matrix

M=zeros(GlobalDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature(2);
 
% cycle for element
for element=1:totalElements       
  % indice : nodal condofectivities for each element
  nodes=elementConnectivityMat(element,:);           
  ndof=length(nodes);
  
  % cycle for Gauss point
  for i=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(i,:);
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [shapeFunctionN,derivativeShapeFunction_xi_eta]=shapeFunction(xi,eta);
    [J,derivativeShapeFunction_x_y]=Jacobian(nodalCoordinateMat(nodes,:),derivativeShapeFunction_xi_eta);
    
    
% mass matrix 
    M(nodes,nodes)=M(nodes,nodes)+shapeFunctionN*shapeFunctionN'*thickness*rho*gaussWeights(i)*det(J);
    
    M(nodes+totalNodes,nodes+totalNodes)=M(nodes+totalNodes,nodes+totalNodes)+shapeFunctionN*shapeFunctionN'*I*rho*gaussWeights(i)*det(J);
    
    M(nodes+2*totalNodes,nodes+2*totalNodes)=M(nodes+2*totalNodes,nodes+2*totalNodes)+shapeFunctionN*shapeFunctionN'*I*rho*gaussWeights(i)*det(J);
    
    M(nodes+3*totalNodes,nodes+3*totalNodes)=M(nodes+3*totalNodes,nodes+3*totalNodes)+shapeFunctionN*shapeFunctionN'*thickness*rho*gaussWeights(i)*det(J);
    
    M(nodes+4*totalNodes,nodes+4*totalNodes)=M(nodes+4*totalNodes,nodes+4*totalNodes)+shapeFunctionN*shapeFunctionN'*thickness*rho*gaussWeights(i)*det(J);
    end 
end   

