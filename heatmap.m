% Box lengths considered
L_x = [32	20	16	10	8];
% Longest wavelength in the box 
for i=1:size(L_x,2)
    m(i) = 1/L_x(i);
end

% top R_m for which mean-field is valid  (R_m << 1)
R_min = 0.4;
% step R_m to consider R_m =  R_min/step_R_min ++ 
step_R_min = 1001;
% set incr=D/D_critical
%%incr = 2.0;
% set D
D = 1.0;


%count plots
k=1;
% for each m(j) calculate the growth rates and plot as heatmaps
figure(1)

for j=1:size(L_x,2)
    clear rm x
    % y-axis is for R_m
    rm(1)= R_min/step_R_min;
    for i=2:step_R_min
        rm(i) = rm(i-1) + R_min/step_R_min;
    end
    
    % x-axis is for x/L, where  x=[0, boxlength_x ], ie x/L=[0,2pi] 
    x(1) = 0.0;
    for i=2:step_R_min; 
        x(i) = (i-1)./(step_R_min-1)*2*pi;
    end
    
    % a mesh because you need data at every possible x-rm combination
    [X,RM] = meshgrid(x,rm);
    % calculate the growth rate, given D-increase
    % p = m^2cos^2(mx)*incr^2/R_min * Rm - 1/Rm * m^2
    %%%% p = m(j)^2.*(cos(X).^2.*incr.^2./R_min.^2.*RM -1./RM); 
    
    % calculate the growth rate, given D itself
    %p = m(j)*( D^2*cos(X).^2.*RM -m(j)./RM); 
    p = D.*m(j).*cos(X).*(RM-m(j)./RM);
    % plot as heatmaps 
    if (k<6)
        clims = [0.0 0.005];
    else
        clims = [0.0 0.005];
    end 
     
    subplot(2*size(L_x,2),1, k);
    k=k+1;
    imagesc(x,rm,p, clims);
    hold on;
    plot(x, rm(end)*cos(x), 'k:','linewidth',1.3);
    plot([0 x(end)],[sqrt(m(j)) sqrt(m(j))],'k-', 'linewidth',1.1);
    hold off;
    % a mesh because
    %[X,Y] = meshgrid(x,y);
    % you need data at every possible x-y combination
    %Z = sqrt(1./X)./Y;
    % that scales the Z-values and plots a heatmap
    %imagesc(x,y,Z)
    set(gca,'xtick',[0,pi/2,pi,3*pi/2, 2*pi])
    set(gca,'xticklabel',{'0';'\pi/2'; '\pi';'3\pi/2';' 2\pi'})
    % choose a colormap of your liking
    colormap autumn(50)
    colorbar
    xlabel('X/L_x')
    ylabel('R_m')
    title(['Growth rate: AR_{box}= {',num2str(L_x(j)),'}'])
    set(gcf,'PaperPositionMode','auto')
    set(gca,'YDir','normal')
%     if (k<7)
%         clims2 = [0.0 0.005];
%     else
%         clims2 = [0.0 0.005];
%     end 
%     
%     subplot(2*size(L_x,2),2, k);
%     k=k+1;
%     contourf(x,rm,p,clims2)
%     hold on;
%     plot(x, rm(end)*cos(x), 'k:');
%     plot(x, sqrt(3.*m(j)),'b:');
%     hold off;
%     set(gca,'xtick',[0,pi/2,pi,3*pi/2, 2*pi])
%     set(gca,'xticklabel',{'0';'\pi/2'; '\pi';'3\pi/2';' 2\pi'})
%     xlabel('X/L_x')
%     ylabel('R_m')
%     title(['Growth rate contour: AR_{box}= {',num2str(L_x(j)),'}'])
%     set(gca,'YDir','normal')
%     colorbar
end 

