

rep1='/store/ASTRO/vs391/kinematic_dynamo/u_abc/results';
rep2='/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/results';

files1=[{[rep1,'/box_8/result_summary']}... 
        {[rep1,'/box_10/result_summary']}...
        {[rep1,'/box_16/result_summary']}...  
        {[rep1,'/box_20/result_summary']}... 
        {[rep1,'/box_32/result_summary']}];
        
legendInfo1 = [     {'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'}];

files2=[{[rep2,'/box_8/result_summary']}... 
        {[rep2,'/box_10/result_summary']}...
        {[rep2,'/box_16/result_summary']}...
        {[rep2,'/box_20/result_summary']}...
        {[rep2,'/box_32/result_summary']}];
       
legendInfo2 = [     {'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'}];
                   
                
color{1} = [0.83,0,0.17]; color{2} = [1,0.67,0.33]; color{3} = [0.33,0.33,1.0];   
color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5]; 
            
full_file=importdata(files1{1,1});
timevar=full_file.data;
n=size(timevar,1);
tblC = cell(100,size(timevar,2));
 for j=1:size(files1,2)
    full_file=importdata(files1{1,j});
    timevar=full_file.data;
    nvar=size(timevar,2);
    nentry=n;
    n=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,2);
    Lx = tblB{:,1};
    Rm = tblB{:,2};
    Gr = tblB{:,3};
    lBx= tblB{:,4};
    lBy= tblB{:,5};
    lBz= tblB{:,6};
    
    % Save data to a complete table
    
    if (j==1)
        for iii=1:nvar
            for ii=1:nentry    
                tblC{ii,iii} = tblB{ii,iii} ; 
            end 
        end
    else
        for iii=1:nvar
            ii1 = 1;
            for ii=(j-1)*(nentry+1): (j-1)*(nentry+1)+n-1 
                tblC{ii,iii} = tblB{ii1,iii} ;
                ii1=ii1+1;
            end 
        end
    end
    figure(1)
        subplot(2,2,1)
            plot(Rm,Gr,'o',...
                'LineStyle', ':',...
                'color',color{j},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',2.5);
            hold on;
            title({'Classic ABC flow'},'fontsize',16, 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo1,'Location','southeast')
            %xlim([0 5]);
            %ylim([0 0.25]);
        subplot(2,2,3)
            plot(Rm,max(lBy,lBz)./(2*pi),'o','LineStyle', ':','color',color{j},'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',2.5);
            hold on;
            title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
            %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
            %legend(legendInfo1,'Location','east')
            ylim([0 32]);
            %line([0 5],[1 1],'LineStyle', ':', 'LineWidth',1.5, 'color', 'k');
%         subplot(3,2,5)
%             plot(lBy./(2*pi),Gr,'o', 'LineStyle', ':','color',color{j},'LineWidth',1.5,...
%                 'MarkerEdgeColor', color{j}, ...
%                 'MarkerFaceColor', color{j}, ...
%                 'MarkerSize',2.5);
%             hold on;
%             %title({'Correlation'},'fontsize',16, 'Interpreter', 'latex');
%             xlabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex');
%             ylabel('$\sigma$','fontsize',16, 'Interpreter', 'latex'  );
%             legend(legendInfo1,'Location','east')
              
 end

  
%clearvars -except files2 legendInfo2
 
full_file=importdata(files2{1,1});
timevar=full_file.data;
n=size(timevar,1);
tblC = cell(100,size(timevar,2));
 for j=1:size(files2,2)
    %color{j} = rand(1,3); 
    full_file=importdata(files2{1,j});
    timevar=full_file.data;
    nvar=size(timevar,2);
    nentry=n;
    n=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,2);
    Lx = tblB{:,1};
    Lx = Lx(2);
    Rm = tblB{:,2}.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2));
    Gr = tblB{:,3};
    lBx= tblB{:,4};
    lBy= tblB{:,5};
    lBz= tblB{:,6};
    
    % Save data to a complete table
    
    if (j==1)
        for iii=1:nvar
            for ii=1:nentry    
                tblC{ii,iii} = tblB{ii,iii} ; 
            end 
        end
    else
        for iii=1:nvar
            ii1 = 1;
            for ii=(j-1)*(nentry+1): (j-1)*(nentry+1)+n-1 
                tblC{ii,iii} = tblB{ii1,iii} ;
                ii1=ii1+1;
            end 
        end
    end
    figure(1)
        subplot(2,2,2)
            plot(Rm,Gr,'o',...
                'LineStyle', ':',...
                'color',color{j},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',2.5);
            hold on;
            title({'Classic modified flow'},'fontsize',16, 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            %ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            %legend(legendInfo2,'Location','east')
            %xlim([0 5]);
            %ylim([0 0.25]);
        subplot(2,2,4)
            plot(Rm,max(lBy,lBz)./(2*pi),'o','LineStyle', ':','color',color{j},'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',2.5);
            hold on;
            title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
            %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            %ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
            %legend(legendInfo2,'Location','east')
            ylim([0 32]);
            %line([0 5],[1 1],'LineStyle', ':', 'LineWidth',1.5, 'color', 'k');
%         subplot(3,2,6)
%             plot(lBy./(2*pi),Gr,'o', 'LineStyle', ':','color',color{j},'LineWidth',1.5,...
%                 'MarkerEdgeColor', color{j}, ...
%                 'MarkerFaceColor', color{j}, ...
%                 'MarkerSize',2.5);
%             hold on;
%             %title({'Correlation'},'fontsize',16, 'Interpreter', 'latex');
%             xlabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex');
%             %ylabel('$\sigma$','fontsize',16, 'Interpreter', 'latex'  );
%             %legend(legendInfo2,'Location','east')
%               
 end
 
 
 for jj=1:5
     clearvars -except jj files1 files2 legendInfo3 legendInfo2
     files3{1,1} = files1{1,jj};    % abc
     files3{1,2} = files2{1,jj};    % iii
     legendInfo3 = [    {'$\vec{u}_{abc}$'}...
                        {'$\vec{u}_{iii}$'}];
     color{1} = [1,0,0];
     color{2} = [0,0,0];
     for j=1:size(files3,2)
        clear timevar Lx Rm Gr lBx lBy lBz
        full_file=importdata(files3{1,j});
        timevar=full_file.data;
        n=size(timevar,1);
        % Create a table with all variables
        %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
        tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,2);
        Lx = tblB{:,1};
        Lx = Lx(2);
        Rm = tblB{:,2};
        Gr = tblB{:,3};
        lBx= tblB{:,4};
        lBy= tblB{:,5};
        lBz= tblB{:,6};
        
        
        
        %fprintf('jj=%d,j=%d\n',jj,j); %, Lx=%d,eta=%f,Rm=%f\n',jj,j,Lx,Rm, Rm*sqrt(0.5 + (1/Lx).^2/3));
        
        
        
%         figure(3)
%         subplot(2,5,jj)
%                 plot(Rm,Gr,'o',...
%                     'LineStyle', ':',...
%                     'color',color{j},...
%                     'LineWidth',1.5,...
%                     'MarkerEdgeColor', color{j}, ...
%                     'MarkerFaceColor', color{j}, ...
%                     'MarkerSize',2.5);
%                 hold on;
%                 title(legendInfo2{1,jj},'fontsize',16, 'Interpreter', 'latex');
%                 xlabel('$1/\eta$','fontsize',16, 'Interpreter', 'latex');
%                 ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
%                 legend(legendInfo3,'Location','southeast','fontsize',16, 'Interpreter', 'latex')
%        subplot(2,5,jj+5)
%                 plot(Rm,max(lBy,lBz)./(2*pi),'o','LineStyle', ':','color',color{j},'LineWidth',1.5,...
%                     'MarkerEdgeColor', color{j}, ...
%                     'MarkerFaceColor', color{j}, ...
%                     'MarkerSize',2.5);
%                 hold on;
%                 %title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
%                 %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
%                 xlabel('$1/\eta$','fontsize',16, 'Interpreter', 'latex');
%                 ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
%                 legend(legendInfo3,'Location','northeast','fontsize',16, 'Interpreter', 'latex')
%                 ylim([0 32]);
%                 
       figure(4)
        subplot(2,5,jj)
                if (mod(j,2)==0)
                    semilogy(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                    %for iii=1:size(Rm)
                    %    fprintf('jj=%d,j=%d,Lx=%d, eta=%f,Rm=%f\n',jj,j,Lx(iii), Rm(iii),Rm(iii)*sqrt(0.5 + (1.0/Lx(iii)).^2/3.) );
                    %end
                else 
                    semilogy(Rm,Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                end
                hold on; set(gca, 'FontSize', 16)
                title(legendInfo2{1,jj},'fontsize',16, 'Interpreter', 'latex');
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                if (jj==1)
                    ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
                    legend(legendInfo3,'Location','southeast','fontsize',16, 'Interpreter', 'latex') 
                end
                ylim([0 0.25]); 
                xlim([0 30]);
       subplot(2,5,jj+5)
                if (mod(j,2)==0)
                    semilogy(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),max(lBz,lBy)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                    ylim([1 32]);
                else 
                    semilogy(Rm,max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                    ylim([1 32]);
                end
                hold on; set(gca, 'FontSize', 16)
                %title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
                %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                if (jj==1)
                ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
                legend(legendInfo3,'Location','northeast','fontsize',16, 'Interpreter', 'latex')
                end
                ylim([1 32]);
                xlim([0 30]);     
        if (jj>2)        
        figure(5)
        subplot(2,3,jj-2)
                if (mod(j,2)==0)
                    semilogy(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                    %for iii=1:size(Rm)
                    %    fprintf('jj=%d,j=%d,Lx=%d, eta=%f,Rm=%f\n',jj,j,Lx(iii), Rm(iii),Rm(iii)*sqrt(0.5 + (1.0/Lx(iii)).^2/3.) );
                    %end
                else 
                   semilogy(Rm,Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                end
                hold on;
                set(gca, 'FontSize', 16)
                title(legendInfo2{1,jj},'fontsize',16, 'Interpreter', 'latex');
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                if (jj==3)
                    ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
                    legend(legendInfo3,'Location','southeast','fontsize',16, 'Interpreter', 'latex');
                end
                xlim([0 2]);
                ylim([-0.1 0.15]);  
       subplot(2,3,jj+1)
                if (mod(j,2)==0)
                    semilogy(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                else 
                    semilogy(Rm,max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',2.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.5);
                end
                hold on;
                %title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
                %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                set(gca, 'FontSize', 16)
                if (jj==3)
                    ylabel('$l(B)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
                    legend(legendInfo3,'Location','northeast','fontsize',16, 'Interpreter', 'latex')
                end
                ylim([0 32]);   
                xlim([0 2]);
        end       
         if (jj==10)
            figure(6)
                 subplot(2,1,1)
                if (mod(j,2)==0)
                    plot(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',2.5);
                    %for iii=1:size(Rm)
                    %    fprintf('jj=%d,j=%d,Lx=%d, eta=%f,Rm=%f\n',jj,j,Lx(iii), Rm(iii),Rm(iii)*sqrt(0.5 + (1.0/Lx(iii)).^2/3.) );
                    %end
                else 
                    plot(Rm,Gr,'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',2.5);
                end
                hold on;
                title(legendInfo2{1,jj},'fontsize',16, 'Interpreter', 'latex');
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                ylabel('$\sigma$ $[T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
                legend(legendInfo3,'Location','southeast','fontsize',16, 'Interpreter', 'latex');
                xlim([0 2]);
       subplot(2,1,2)
                if (mod(j,2)==0)
                    plot(Rm.*(1.5+1./Lx.^2)/sqrt(1.5*(1+2*1./Lx.^4+9*1./Lx.^2)),max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',2.5);
                else 
                    plot(Rm,max(lBy,lBz)./(2*pi),'o',...
                        'LineStyle', ':',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{j}, ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',2.5);
                end
                hold on;
                %title({'Length scale of the mean-field'},'fontsize',16, 'Interpreter', 'latex');
                %plot(Rm,lBz./(2*pi),'LineStyle', '-','color',color{j},'LineWidth',1.5);
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                ylabel('$l(B_y)/2\pi$ $[L]$','fontsize',16, 'Interpreter', 'latex'  );
                %legend(legendInfo3,'Location','northeast','fontsize',16, 'Interpreter', 'latex')
                ylim([0 32]);   
                %xlim([0 2]);
         end
               
     end         
 end  
          
            