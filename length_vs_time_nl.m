close all;
clear all;

box=[10 32];

color{1}=[0.83 0.17 0.17];
color{2}=[0 0 0];

for jj=1:size(box,2)
    clear emf evf ratio
    rep=['/store/ASTRO/vs391/nonlinear_dynamo/u_iii/box_', num2str(box(jj))];
    file_info=[rep, '/info'];
    full_file=importdata(file_info);
    info=full_file;
    tblA = table(info(:,1),info(:,2), info(:,3), info(:,4));
    % Sort the rows of the table based on Re, then Rm
    tblB = sortrows(tblA,[4,3]); 
    run = tblB{1:end,1}; 
    cas = tblB{1:end,2}; 
    Rmm  = tblB{1:end,3}; 
    Ree  = tblB{1:end,4}; 
    
    countb=0;
    countu=1;
    clear lebx leby lebz leux leuy leuz;
    for j=1:size(Ree,1)
        file_time_b=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar2'];
        full_file_b=importdata(file_time_b);
        file_time_u=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar3'];
        full_file_u=importdata(file_time_u);
        timevar_b=full_file_b;
        timevar_u=full_file_u;
        
        tb  = timevar_b(:,1);
        lbx = timevar_b(:,5);
        lby = timevar_b(:,6);
        lbz = timevar_b(:,7);
        tu  = timevar_u(:,1);
        lux = timevar_u(:,5);
        luy = timevar_u(:,6);
        luz = timevar_u(:,7);
    %%%%%%%%%%
    %% PLOT %%
    %%%%%%%%%%
    kk=4;  % number of columns 
    % num is number of rows
    if mod(size(Ree,1),kk) > 0
        num = fix(size(Ree,1)/kk) + 1; 
    else
        num = fix(size(Ree,1)/kk);
    end
    hFig = figure(jj);
    set(hFig, 'Position', [100, 50, 1100, 500]);
    fprintf( 'j=%d, cb=%d\n',j, countb)
    subplot(2*num,kk,j+countb*kk)
        p1= plot(tb,lbx/6.28, ...
            'LineStyle', '-',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_x^B');
        ylim([0 box(jj)])
        %if countb==0
            title(['R_e=',num2str(Ree(j)),', R_m=',num2str(Rmm(j))] )
        %end
        hold on;
        p2= plot(tb,lby/6.28, ...
            'LineStyle', '-.',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_y^B');
        p3= plot(tb,lbz/6.28, ...
            'LineStyle', ':',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_z^B');
        if mod(j,kk)==0
            countb=countb+1;
        end
        if j==size(Ree,1)
            legend([p1 p2 p3],'Location','southwest');
        end
     subplot(2*num,kk,j+countu*kk)   
        p4= plot(tu,lux/6.28, ...
            'LineStyle', '-',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_x^u'); 
        ylim([0 box(jj)])
        if countu==num*2-1
            xlabel('Time')
        end
        hold on;
        p5= plot(tu,luy/6.28, ...
            'LineStyle', '-.',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_y^u');
        
        p6= plot(tu,luz/6.28, ...
            'LineStyle', ':',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',1.5,...
            'DisplayName','L_z^u');
        
        if mod(j,kk)==0
            countu=countu+1;
        end
        if j==size(Ree,1)
            legend([p4 p5 p6],'Location','southwest');
        end
        
        lebx(j)=lbx(end);
        leby(j)=lby(end);
        lebz(j)=lbz(end);
        leux(j)=lux(end);
        leuy(j)=luy(end);
        leuz(j)=luz(end);
    end
    
    hFig1 = figure(jj+size(box,2));
    set(hFig1, 'Position', [100, 50, 600, 500]);
        p7=plot(Rmm./Ree, lebx, 'o',...
            'LineStyle', '-',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_x^B'); 
        xlabel('Pm')
        ylabel('[L]')
        hold on;
        p8=plot(Rmm./Ree, leby, 'o',...
            'LineStyle', '-.',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_y^B'); 
        p9=plot(Rmm./Ree, lebz, 'o',...
            'LineStyle', ':',...
            'color',color{1},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_z^B');      
        p10=plot(Rmm./Ree, leux, 'o',...
            'LineStyle', '-',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_x^u'); 
        hold on;
        p11=plot(Rmm./Ree, leuy, 'o',...
            'LineStyle', '-.',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_y^u'); 
        p12=plot(Rmm./Ree, leuz, 'o',...
            'LineStyle', ':',...
            'color',color{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',3.5,...
            'DisplayName','L_z^u');    
        legend([p7 p8 p9 p10 p11 p12],'Location','west');
end