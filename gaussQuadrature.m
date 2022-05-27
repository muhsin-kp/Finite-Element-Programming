function [gaussWeights,gaussPoints]=gaussQuadrature(Npoint)

    if Npoint==2
            gaussPoints=1/sqrt(3)*[ -1 -1
                                        1 -1
                                        1  1
                                       -1  1];
            gaussWeights=[ 1
                    1
                    1
                    1]; 
    elseif Npoint==1
            gaussPoints=[0 0];
            gaussWeights=4;
    end
end