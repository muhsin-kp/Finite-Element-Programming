

function [K]=calculateStiffnessMatrix(GlobalDof,totalElements,...
    elementConnectivityMat,totalNodes,nodalCoordinateMat,AMatrix,...
    BMatrix,DMatrix,SMatrix)

    % Initiate Global Stiffness matrix
    
    K=zeros(GlobalDof);
    
    [gaussWeights,gaussLocations]=gaussQuadrature(2);
     
    % For Bending Stifness matrix- loop through element
    for element=1:totalElements  
        
      nodes=elementConnectivityMat(element,:);           
      nodalDOF=[ nodes nodes+totalNodes nodes+2*totalNodes ...
                nodes+3*totalNodes nodes+4*totalNodes]; 
      n=length(nodes);
      
      for i=1:size(gaussWeights,1)                      
        gaussPoint=gaussLocations(i,:);                                                     
        xi=gaussPoint(1);
        eta=gaussPoint(2);
    
    % shape functions and derivatives and jacobian
    
        [shapeFunctionN,derivativeShapeFunction_xi_eta]=shapeFunction(xi,eta);
    
        [J,derivativeShapeFunction_x_y]=Jacobian(nodalCoordinateMat(nodes,:),derivativeShapeFunction_xi_eta);
        
    % B matrix bending
        Bbending=zeros(3,5*n);
        Bbending(1,n+1:2*n)        = derivativeShapeFunction_x_y(1,:);
        Bbending(2,2*n+1:3*n)      = derivativeShapeFunction_x_y(2,:);
        Bbending(3,n+1:2*n)        = derivativeShapeFunction_x_y(2,:);
        Bbending(3,2*n+1:3*n)      = derivativeShapeFunction_x_y(1,:);
    % B matrix membrane
        Bmembrane=zeros(3,5*n);
        Bmembrane(1,3*n+1:4*n)      = derivativeShapeFunction_x_y(1,:);
        Bmembrane(2,4*n+1:5*n)      = derivativeShapeFunction_x_y(2,:);
        Bmembrane(3,3*n+1:4*n)      = derivativeShapeFunction_x_y(2,:);
        Bmembrane(3,4*n+1:5*n)      = derivativeShapeFunction_x_y(1,:);
        
    
    % bending-bending stiffness matrix
        K(nodalDOF,nodalDOF)=K(nodalDOF,nodalDOF)+...
                           Bbending'*DMatrix*Bbending*gaussWeights(i)*det(J);
    % membrane-membrane stiffness matrix                 
        K(nodalDOF,nodalDOF)=K(nodalDOF,nodalDOF)+...
                           Bmembrane'*AMatrix*Bmembrane*gaussWeights(i)*det(J);
    % membrane-bending stiffness matrix                 
        K(nodalDOF,nodalDOF)=K(nodalDOF,nodalDOF)+...
                           Bmembrane'*BMatrix*Bbending*gaussWeights(i)*det(J);
    % bending-membrane stiffness matrix                
        K(nodalDOF,nodalDOF)=K(nodalDOF,nodalDOF)+...
                           Bbending'*BMatrix*Bmembrane*gaussWeights(i)*det(J);
    
      end
    end   
    
    % shear stiffness matrix
    
    % Gauss quadrature for shear part
    [gaussWeights,gaussLocations]=gaussQuadrature(1);
     
    % loop through element
    for element=1:totalElements       
      
      nodes=elementConnectivityMat(element,:);           
      nodalDOF=[ nodes nodes+totalNodes nodes+2*totalNodes ...
                nodes+3*totalNodes nodes+4*totalNodes]; 
      n=length(nodes);
    
      for i=1:size(gaussWeights,1)                      
        gaussPoint=gaussLocations(i,:);                                                     
        xi=gaussPoint(1);
        eta=gaussPoint(2);
    
    % shape functions and derivatives
        
        [shapeFunctionN,derivativeShapeFunction_xi_eta]=shapeFunction(xi,eta);
    
    % Jacobian matrix, inverse of Jacobian, 
    
        [J,derivativeShapeFunction_x_y]=Jacobian(nodalCoordinateMat(nodes,:),derivativeShapeFunction_xi_eta);
        
    % Bshear Matrix
    
        Bshear=zeros(2,5*n);
        Bshear(1,1:n)       = derivativeShapeFunction_x_y(1,:);
        Bshear(2,1:n)       = derivativeShapeFunction_x_y(2,:);
        Bshear(1,n+1:2*n)  = shapeFunctionN;
        Bshear(2,2*n+1:3*n)= shapeFunctionN;
    
    % stiffness matrix shear
        K(nodalDOF,nodalDOF)=K(nodalDOF,nodalDOF)+...
            Bshear'*SMatrix  *Bshear*gaussWeights(i)*det(J);  
      end
end 
