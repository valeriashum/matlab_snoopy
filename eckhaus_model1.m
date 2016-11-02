clear all;
close all;

eta=[0.01 1.0]
for ii = 1:2
    for i=1:2
        %eta=i*10.0;
        x=0:0.01:1.0; %k
        y=0.1:0.01:1.1;%q
        z=1:0.1:11; %RO
        [X,Y,Z] = meshgrid(x,y,z);
        if ii==1 %(model1)
            v = -eta * X.^2 + X.*(1-2*eta.*Y - 4*Z.^2 + 3*Z.^4) + Y.*(-eta*Y + (1-3*Z.^2).^2);
        % f= -eta k^2 + k (1 - 2 eta q - 4 R0^2 + 3 R0^4) + 
        %    q (-eta q + (1 - 3 R0^2)^2);
        else %(model2)
            v = -eta*(X + Y).^2 + (X + Y).*(3*Y.*Z + 1).^2;
        end

        figure(ii)
        subplot(1,2,i)
            scatter3(X(:), Y(:), Z(:), 72, v(:), 'filled');
            hold on
            view(-15, 35);
            colorbar;
            xlabel('k')
            ylabel('q')
            zlabel('R_0')
        if ii==1
            title('Model 1: maximum growth rate')
        else
            title('Model 2: maximum growth rate')
        end
    end
end