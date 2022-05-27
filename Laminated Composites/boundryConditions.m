
function [unConstrainedDOFs,constrainedNodesinW]=...
    boundryConditions(BC,GlobalDof,nodalCoordinateMat,totalNodes)

    xx=nodalCoordinateMat(:,1);
    yy=nodalCoordinateMat(:,2);
    switch BC
        case 's-s-s-s'
            constrainedNodesinW =find(xx==max(nodalCoordinateMat(:,1))|...
                             xx==min(nodalCoordinateMat(:,1))|...
                             yy==min(nodalCoordinateMat(:,2))|...
                             yy==max(nodalCoordinateMat(:,2)));
           
            constrainedNodesinThetaX =find(yy==max(nodalCoordinateMat(:,2))|...
                             yy==min(nodalCoordinateMat(:,2)));

            constrainedNodesinThetaY =find(xx==max(nodalCoordinateMat(:,1))|...
                             xx==min(nodalCoordinateMat(:,1)));

            constrainedNodesinU =find(xx==min(nodalCoordinateMat(:,1)));

            constrainedNodesinV =find(yy==min(nodalCoordinateMat(:,2)));

        case 'c-c-c-c'
    
            constrainedNodesinW =find(xx==max(nodalCoordinateMat(:,1))|...
                             xx==min(nodalCoordinateMat(:,1))|...
                             yy==min(nodalCoordinateMat(:,2))|...
                             yy==max(nodalCoordinateMat(:,2)));

            constrainedNodesinThetaX =constrainedNodesinW;

            constrainedNodesinThetaY =constrainedNodesinThetaX;

            constrainedNodesinU =constrainedNodesinThetaX;

            constrainedNodesinV =constrainedNodesinThetaX;
    
    end
    
    constrainedDOFs=[constrainedNodesinW;constrainedNodesinThetaX+totalNodes;...
          constrainedNodesinThetaY+2*totalNodes;...
          constrainedNodesinU+3*totalNodes;constrainedNodesinV+4*totalNodes];

    unConstrainedDOFs=setdiff([1:GlobalDof]',[constrainedDOFs]);