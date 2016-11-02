close all;
clear all;

rep={'/store/ASTRO/vs391/kinematic_dynamo/u_abc/results', '/store/ASTRO/vs391/kinematic_dynamo/u_iii/results'};
box=[8 10 16 20 32];
color{1} = [0.83 0.17 0.17];
color{2} = [0 0 0];
color{3} = [0.33,0.33,1.0];   
color{4} = [0,0.5,0]; 
color{5} = [0.5,0,0.5]; 

for ii = 1:size(rep,2)
    for jj = 1:size(box,2)
        clear lBx lBy lBz Gr Rm Lx
        legendInfo1{jj} = ['Box=[', num2str(box(jj)), ',2,1]']; 
        file_time=[rep{ii},'/box_', num2str(box(jj)), '/result_summary'];
        full_file=importdata(file_time);
        timevar=full_file.data;
        tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,2);
        Lx = tblB{:,1};
        if ii == 1
            Rm = tblB{:,2};
        else
            Rm = tblB{:,2}.*(1.5+1./Lx.^2)/sqrt(1.5*(3+2*1./Lx.^4+9*1./Lx.^2));
        end
        Gr = tblB{:,3};
        lBx= tblB{:,4};
        lBy= tblB{:,5};
        lBz= tblB{:,6};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hfig = figure(1);
        set(hFig, 'Position', [100, 50, 1100, 900]); 
        set(gca, 'FontSize', 14)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,2,ii)
            plot(Rm,Gr,'o',...
                'LineStyle', ':',...
                'color',color{jj},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{jj}, ...
                'MarkerFaceColor', color{jj}, ...
                'MarkerSize',2.5);
            hold on;
            if ii ==1
                title({'u_{ABC}'},'fontsize',14, 'Interpreter', 'latex');
            else 
                title({'u_{m}'},'fontsize',14, 'Interpreter', 'latex');
            end
            if jj == 1
                ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex'); 
                legend(legendInfo1,'Location','southeast')
            end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,2,ii + 2)
            plot(Rm,max(lBy,lBz)./(2*pi),'o','LineStyle', ':','color',color{j},'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',2.5);
            hold on;
            if ii == 2
                xlabel('$R_m$','fontsize',14, 'Interpreter', 'latex');
            end
            if jj == 1
            ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
            end
            ylim([0 32]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        hfig = figure(4);
        set(hFig, 'Position', [100, 50, 1100, 900]); 
        set(gca, 'FontSize', 14)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,5,jj)
            semilogy(Rm,Gr,'o',...
                'LineStyle', ':',...
                'color',color{ii},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{ii}, ...
                'MarkerFaceColor', color{ii}, ...
                'MarkerSize',2.5);
            hold on; 
            title(['Box=[', num2str(box(jj)), ',2,1]']);
            
                if (jj==1)
                    ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
                    legend(legendInfo3,'Location','southeast','fontsize',16, 'Interpreter', 'latex') 
                end
                ylim([0 0.25]); 
                xlim([0 30]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        subplot(2,5,jj+5)
            semilogy(Rm,max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',2.5);
                hold on; 
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                if (jj==1)
                ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
                legend(legendInfo3,'Location','northeast','fontsize',16, 'Interpreter', 'latex')
                end
                ylim([1 32]);
                xlim([0 30]);     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hfig = figure(6);
        set(hFig, 'Position', [100, 50, 1100, 900]); 
        set(gca, 'FontSize', 14) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,1,1)
            plot(Rm,Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{jj},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{jj}, ...
                        'MarkerFaceColor', color{jj}, ...
                        'MarkerSize',2.5);
                hold on;
                title(['Box=[', num2str(box(jj)), ',2,1]']);
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                ylabel('$\sigma$ $[T^{-1}_{turnover}]$', 'Interpreter', 'latex')
                xlim([0 2]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,1,2)
            plot(Rm,max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{jj},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{jj}, ...
                        'MarkerFaceColor', color{jj}, ...
                        'MarkerSize',2.5);
                hold on;
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
                ylim([0 32]);   
                %xlim([0 2]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
     end         
 end  
          
            