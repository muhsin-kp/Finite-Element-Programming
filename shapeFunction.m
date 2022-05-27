function [shapeFunctions,derivativeShapeFunction_si_eta]=shapeFunction(si,eta)

    % Shape function for quadrilateral isoparamteric element

    shapeFunctions=[ (1-si)*(1-eta)/4
            (1+si)*(1-eta)/4
            (1+si)*(1+eta)/4
            (1-si)*(1+eta)/4];
    
    % Derivative of Shape function for quadrilateral isoparamteric element
    % w.r.t si and eta

    derivativeShapeFunction_si_eta=[-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4
                        -(1-si)/4,-(1+si)/4,(1+si)/4,(1-si)/4];

end