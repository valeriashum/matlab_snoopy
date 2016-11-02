close all

rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

box=[8 10 16 20 32];
kkk=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

files1=[{[rep1,'/b_meas_8_copy.dat']}...
        {[rep1,'/b_meas_10_copy.dat']}...
        {[rep1,'/b_meas_16_copy.dat']}...
        {[rep1,'/b_meas_20_copy.dat']}...
        {[rep1,'/b_meas_32_copy.dat']} ];   

files2=[{[rep2,'/b_meas_8_copy.dat']}...
        {[rep2,'/b_meas_10_copy.dat']}...
        {[rep2,'/b_meas_16_copy.dat']}...
        {[rep2,'/b_meas_20_copy.dat']}...
        {[rep2,'/b_meas_32_copy.dat']} ];    
    
legendInfo1 = [     {'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'}];
               
color{1} = [0.83,0,0.17]; color{2} = [0,0,0];               
for ind=1:10
    color{ind}=rand(1,3);
end

for FLOW=1:2
    if FLOW ==1
        files = files1; 
    else
        files = files2;
    end

% FIGURE 1
for j=1:size(files,2)

    full_file=importdata(files{1,j});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3), ...
                timevar(:,4),timevar(:,5),timevar(:,6), ...
                timevar(:,7),timevar(:,8),timevar(:,9), ...
                timevar(:,10),timevar(:,11),timevar(:,12));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,1); 
    if FLOW==1
        Lx = box(j);
        Rm   = tblB{:,1}*(1.5+1./Lx.^2)/sqrt(1.5*(3+2*1./Lx.^4+9*1./Lx.^2));
    else 
        Rm   = tblB{:,1};
    end
    LB   = tblB{:,2};
    for ind=1:10
        clear Ratio;
        Ratio = tblB{:,2+ind};
        figure(1)
        set(gca, 'FontSize', 14)
            subplot(1,6,j)
                loglog(Rm,Ratio,'o',...
                        'LineStyle', 'none',... 
                        'color',color{ind},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ind}, ...
                        'MarkerFaceColor', color{ind}, ...
                        'MarkerSize',4.5);     
             hold on;
   
%     if (j==5)
%     legend({'$k_x=0.1$','$k_x=0.2$',...
%            '$k_x=0.3$','$k_x=0.4$',...
%            '$k_x=0.5$','$k_x=0.6$',...
%            '$k_x=0.7$','$k_x=0.8$',...
%            '$k_x=0.9$','$k_x=1.0$'},...
%            'Interpreter','latex');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig = figure(10 + j);
        set(hFig, 'Position', [10, 50, 500, 400]);
        if ind == 5
            loglog(Rm,Ratio,'o',...
                        'LineStyle', '-',... 
                        'color',color{ind},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ind}, ...
                        'MarkerFaceColor', color{ind}, ...
                        'MarkerSize',4.5)
            title(legendInfo1{1,j},'fontsize',16);
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$ E_m^{large}/E_m^{small}$','fontsize',16, 'Interpreter', 'latex');
hold on;
            plot([0.1 40],[1 1],'k:')
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
            subplot(1,6,j)
            p = polyfit(log(Rm(Rm<1)),log(Ratio(Rm<1)),1); 
            m1= p(1);
            b1 = exp(p(2));
            loglog(Rm, b1*Rm.^m1, ...
                    'LineStyle', ':',...
                    'color','k',...
                    'LineWidth',1.5);
        Power_laws(j,ind)=m1;
        
        % Calculate interception with ratio = 1
        for jj=1:size(Ratio,1)
            if (Ratio(jj) < 1) 
                y1=Ratio(jj);
                x1=Rm(jj);
                if (jj~=size(Ratio,1))
                    y2=Ratio(jj+1);
                    x2=Rm(jj+1);   
                    if (x1~=0 && y1~=0 && x2~=0 && y2~=0)
                        inters(ind) = (1-y2)*(x2-x1)/(y2-y1) + x2;
                    else 
                        inters(ind)=0;
                    end
                else 
                    inters(ind)=Rm(jj);    
                end
            else
                x1=0;
                x2=0;
                y1=0;
                y2=0;
            end 
        end
        fprintf('x1=%f,x2=%f,y1=%f,y2=%f\n',x1,x2,y1,y2);
        
        
        figure(2)
            subplot(1,5,j)
            set(gca, 'FontSize', 14)
            plot(kkk(inters>0),transpose(inters(inters>0)),'o',...
                        'LineStyle', 'none',... 
                        'color',color{ind},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ind}, ...
                        'MarkerFaceColor', color{ind}, ...
                        'MarkerSize',4.5);
            hold on;
    %ylim([-7.5 -2.5])
    xlabel('$k^{lim}$','Interpreter','latex','fontsize',16)
    %title('Power Law $R_m^n$','Interpreter','latex','fontsize',16)
    ylabel('$R_m^{tr}$','Interpreter','latex','fontsize',16)  
    end
    figure(1)
     title(legendInfo1{1,j},'fontsize',16, 'Interpreter', 'latex');
     xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
     %.*(1.5+1./Lx.^2)/sqrt(1.5*(3+2*1./Lx.^4+9*1./Lx.^2))
     if (j==1)
        ylabel('$ E_m^{large}/E_m^{small}$','fontsize',16, 'Interpreter', 'latex');
     end 
            plot([0.1 100],[1 1],'k:')
            plot([1 1],[0.00000001 1e5],'k:')
            plot([2 2],[0.00000001 1e5],'k:')
            ylim([1e-4 1e5])
            %xlim([0.1 2])
            
   
end

% FIGURE 1
figure(1)

for ind=1:10
    subplot(1,6,6)
    set(gca, 'FontSize', 14)
    plot(box,Power_laws(:,ind),'o',...
                        'LineStyle', ':',... 
                        'color',color{ind},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ind}, ...
                        'MarkerFaceColor', color{ind}, ...
                        'MarkerSize',4.5);
    hold on;
    ylim([-7.5 -2.5])
    xlabel('$L_x$','Interpreter','latex','fontsize',16)
    title('Power Law $R_m^n$','Interpreter','latex','fontsize',16)
    ylabel('$n$','Interpreter','latex','fontsize',16)
end 
end

legend({   '$0.1$','$0.2$',...
           '$0.3$','$0.4$',...
           '$0.5$','$0.6$',...
           '$0.7$','$0.8$',...
           '$0.9$','$1.0$'},...
           'Interpreter','latex',...
           'Orientation','horizontal',...
           'Location','north');





