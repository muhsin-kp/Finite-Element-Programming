function [] = plotModeShape(eigValues,eigVectors,thicknessZ,reqModes,totalNodes,LengthL,LengthW,ElementsX,ElementsY,constrainedNodesinW,materialProperties)
    
    %Assign Material Properties

        EY=materialProperties(2);nyuXY=materialProperties(6);nyuYX=materialProperties(7);
        factor=1-nyuXY*nyuYX;rho=materialProperties(8);
   
    % natural frequency=eigen value, omegabar=non dimensional natural
    % frequency

    D0=EY*thicknessZ^3/12/(1-nyuXY*nyuYX);
    eigValues = diag(sqrt(eigValues));%*LengthW*LengthW/pi/pi*sqrt(rho*thicknessZ/D0));
    omegabar = eigValues*LengthW*LengthW/pi/pi*sqrt(rho*thicknessZ/D0);
    [eigValues,ind] = sort(eigValues); eigValues =eigValues/(2*pi);
    ind = ind(1:reqModes); 
    omegabar=sort(omegabar);
    reqV = eigVectors(:,ind);
    unConstrainedDOFW=setdiff([1:totalNodes]',[constrainedNodesinW]);
    NumberunConstrainedDOF=size(unConstrainedDOFW);
    
    modes(1:totalNodes,1:reqModes)=0;
    for i=1:reqModes
        modes(unConstrainedDOFW,i)=reqV(1:NumberunConstrainedDOF,i);
    end
 
    disp(eigValues(1:reqModes));

    x=linspace(-LengthL,LengthL,ElementsX+1);
    y=linspace(-LengthW,LengthW,ElementsY+1);   
    figure("Name","Mode Shapes","Color",'w')
        for kk=1:reqModes
            subplot(2,reqModes/2,kk)
            MODE=modes(:,kk);
            MODE=reshape(MODE,ElementsX+1,ElementsY+1);
            contourf(x,y,MODE);
            colormap(white);
            %title({strcat(strcat(strcat('MODE :',num2str(kk)),' NF: '),num2str(eigValues(kk)))});
            title("MODE "+num2str(kk)+' NF: '+num2str(eigValues(kk))+" NDF: "+omegabar(kk));
            hold on;
            surf(x,y,MODE/max(abs(MODE(:))));
            colormap(jet);
            hold on;
            shading interp
            colorbar('east')
            xlabel('clamped edge / chord')
            ylabel('span length')
            zlabel('deflection w')
            view(45,45)
            axis equal
            daspect([1 1 3])
            rotate3d
            grid off
            axis off
        end

     figure("Name","Mesh","Color",'w')
        X = x;
        Y= y;
        [X,Y] = meshgrid(X,Y)
        plot(X,Y,'k-');
        hold on 
        plot(Y,X,'k-');
        X = reshape(X',[],1); Y = reshape(Y',[],1);
        text(X, Y, cellstr(num2str((1:(ElementsX+1)*(ElementsY+1))')),"FontSize",8,"VerticalAlignment","top");

        x=linspace(-LengthL,LengthL,ElementsX+1);
        y=linspace(-LengthW,LengthW,ElementsY+1);
        x=conv(x,[0.5,0.5],'valid');
        y=conv(y,[0.5,0.5],'valid');
        X1 = x;
        Y1= y;
        [X1,Y1] = meshgrid(X1,Y1)
        X1 = reshape(X1',[],1); Y1 = reshape(Y1',[],1);
        text(X1, Y1, cellstr(num2str((1:(ElementsX*ElementsY))')),"FontSize",10,"Color",'w','HorizontalAlignment','center');
        daspect([1 1 1])
        set(gca,'XTick',[],'YTick',[],'color','#4DBEEE')