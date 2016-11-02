close all;
clear all;

rep1=[{'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/box_', '/store/ASTRO/vs391/nonlinear_dynamo/u_iii/box_'}];
rep2=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc'     , '/store/ASTRO/vs391/kinematic_dynamo/u_iii'}];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SET ALL_PLOTS=1 to get all plots    %%
        % SET ALL_PLOTS=0 to hide state plots %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
ALL_PLOTS = 1;              %
KINEMATIC_CF_PLOTS = 0;     %
LENGTH_3D_PLOTS = 0;        %
LENGTH_3D_RE1_PLOTS = 0;    %
ENERGY_BOUND_PLOTS = 0;     %
RE1_PLOTS = 0;              %
RATIO_PLOT = 0;             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box=[10 32];
color{1}=[0.83 0.17 0.17];
color{2}=[0 0 0];

count_plots = 0;
count_to_4= 0;
for ii=1:1%size(rep1,2)
    clear tblA tblB run cas Ree emf evf lbf lvf
    
    for jj=1:size(box,2)
        count_ree_1 = 0;
        count_ree_10 = 0;
        count_to_4 = count_to_4 + 1;
        count_ree1 = 0;
        clear emf evf ratio1 ratio2 lbf  lvf ratio_l1 ratio_l2 Ree1 Ree2 Rmm1 Rmm2 pre1 pre2 Rmmkr lbxkr lbykr lbzkr ratio1 ratio2 ratio_l1 ratio_l2   
        rep=[rep1{ii}, num2str(box(jj))];
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
        count1=0;
        count2=0;
        for j=1:size(Ree,1)
        if cas(j) ~= 0
            count_plots = count_plots + 1;
            clear lbx lby lbz lux luy luz timeb timeu timeb2 timeu2 lb2 lu2
            %file_time=[rep,'/run_', num2str(run(j)), '/previous_run_', num2str(cas(j)), '/timevar'];
            file_time=[rep,'/run_', num2str(run(j)), '/for_analysis/timevar'];
            full_file=importdata(file_time);
            timevar=full_file.data;
            t  =timevar(:,1);
            Rm =timevar(:,2);
            Re =timevar(:,3);
            em =timevar(:,4);
            ev =timevar(:,5);
            hv =timevar(:,6);
            w2 =timevar(:,7);
            hc =timevar(:,8);
            hm =timevar(:,9);
            
%             if ii == 2 
%                 m = 1.0/box(jj);
%                 Rmm = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Rmm;
%                 Ree = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Ree;
%             end     

            if run(j) < 100 
                file_time_b=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar2'];
                file_time_u=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar3'];

                file_time_b2=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar4'];
                file_time_u2=[rep,'/run_', num2str(run(j)), '/for_analysis/result1/timevar5'];

                color{1}=[0.83 0.17 0.17];
                color{2}=[0 0 0];
            else
                file_time_b=[rep,'/run_', num2str(run(j)), '/for_analysis/graphics/timevar2'];
                file_time_u=[rep,'/run_', num2str(run(j)), '/for_analysis/graphics/timevar3'];

                file_time_b2=[rep,'/run_', num2str(run(j)), '/for_analysis/graphics/timevar4'];
                file_time_u2=[rep,'/run_', num2str(run(j)), '/for_analysis/graphics/timevar5'];

                color{1}=[1.0 0.8 0.6];
                color{2}=[0.4 0.4 0.4];
            end

            full_file_b=importdata(file_time_b);
            full_file_u=importdata(file_time_u);
            full_file_b2=importdata(file_time_b2);
            full_file_u2=importdata(file_time_u2);
            timevar_b=full_file_b;
            timevar_u=full_file_u;
            timevar_b2=full_file_b2;
            timevar_u2=full_file_u2;

            timeb   = timevar_b(:,1);
            lbx = timevar_b(:,5);
            lby = timevar_b(:,6);
            lbz = timevar_b(:,7);
            timeu   = timevar_u(:,1);
            lux = timevar_u(:,5);
            luy = timevar_u(:,6);
            luz = timevar_u(:,7);

            timeb2   = timevar_b2(:,1);
            lb2 = timevar_b2(:,2);
            b01 = timevar_b2(:,3);
            b02 = timevar_b2(:,4);
            b03 = timevar_b2(:,5);
            b04 = timevar_b2(:,6);
            b05 = timevar_b2(:,7);
            b06 = timevar_b2(:,8);
            b07 = timevar_b2(:,9);
            b08 = timevar_b2(:,10);
            b09 = timevar_b2(:,11);
            b10 = timevar_b2(:,12);

            timeu2   = timevar_u2(:,1);
            lu2 = timevar_u2(:,2);
            u01 = timevar_u2(:,3);
            u02 = timevar_u2(:,4);
            u03 = timevar_u2(:,5);
            u04 = timevar_u2(:,6);
            u05 = timevar_u2(:,7);
            u06 = timevar_u2(:,8);
            u07 = timevar_u2(:,9);
            u08 = timevar_u2(:,10);
            u09 = timevar_u2(:,11);
            u10 = timevar_u2(:,12);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            % if for some reason the lengthscale files arent of 
            % equal size, match them 
            time_l = timeb;
            time_l2 = timeb2;
            %fprintf('j=%d,BEFORE: tb=%d,tu=%d, tb2=%d,tu2=%d\n',j,size(timeb,1), size(timeu,1),size(timeb2,1),size(timeu2,1))
            if size(timeu,1) > size(timeb,1)
                luxx = lux(1:size(timeb,1));      luyy = luy(1:size(timeb,1));      luzz = luz(1:size(timeb,1));
                clear lux luy luz time_l
                lux = luxx;                     luy = luyy;                     luz = luzz;
                time_l = timeb;
                clear luxx luyy luzz
            elseif size(timeu,1) < size(timeb,1)
                lbxx = lbx(1:size(timeu,1));      lbyy = lby(1:size(timeu,1));      lbzz = lbz(1:size(timeu,1));
                clear lbx lby lbz time_l
                lbx = lbxx;                     lby = lbyy;                     lbz = lbzz;
                time_l = timeu;
                clear lbxx lbyy lbzz
            end

             if size(timeu2,1) > size(timeb2,1)
                lu22 = lu2(1:size(timeb2,1));
                clear lu2 time_l2
                time_l2 = timeb2;
                lu2 = lu22;
                clear lu22
             elseif size(timeu2,1) < size(timeb2,1)
                lb22 = lb2(1:size(timeu2,1)); 
                clear lb2 time_l2
                time_l2 = timeu2;
                lb2 = lb22;
                clear lb22
             end
             %fprintf('j=%d, AFTER: tb=%d,tu=%d, tb2=%d,tu2=%d\n',j,size(timeb,1), size(timeu,1),size(timeb2,1),size(timeu2,1))
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % PLOTS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            if ALL_PLOTS == 1
                hFig = figure(count_plots);
                set(hFig, 'Position', [100, 50, 1100, 900]); 
                    subplot(4,5,1:4)
                        p1= semilogy(t,em, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','E_m'); 
                        title(['run=', num2str(run(j)),'     R_e=', num2str(Ree(j)),'        R_m=', num2str(Rmm(j)),'      [', num2str(box(jj)), ',2,1]']);
                        xlim([0 t(end)])
                        ylim([min(em(end),ev(end))*0.8 max(em(end),ev(end))*1.2])
                        hold on;
                        p2= semilogy(t,ev, ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','E_k');
                        legend([p1 p2],'Location','west', 'Orientation','Horizontal')
                        plot([t(1) t(end)],[em(end) em(end)],'k:');
                        plot([t(1) t(end)],[ev(end) ev(end)],'k:');
                   subplot(4,5,5)    
                        plot(t,em./ev, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5);
                            title('E_m/E_k');
                            xlim([0 t(end)])

                    subplot(4,5,6:9)
                        p3= plot(time_l,lby, ...
                            'LineStyle', '--',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^B_y'); 
                        xlim([0 t(end)])
                        hold on;
                        p4= plot(time_l,lbz, ...
                            'LineStyle', '-.',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^B_z'); 
                        p41= plot(time_l2,lb2, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^B'); 

                        p5= plot(time_l,luy, ...
                            'LineStyle', '--',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^u_y'); 
                        p6= plot(time_l,luz, ...
                            'LineStyle', '-.',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^u_y');              

                         p7= plot(time_l2,lu2, ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^u'); 
                        legend([p5 p6 p7 p3 p4 p41],'Location','north', 'Orientation','Horizontal')

                        plot([timeb(1) timeb(end)],[lby(end) lby(end)],'k:');
                        plot([timeb(1) timeb(end)],[lbz(end) lbz(end)],'k:');
                        plot([timeb2(1) timeb2(end)],[lb2(end) lb2(end)],'k:');
                        plot([timeu2(1) timeu2(end)],[lu2(end) lu2(end)],'k:');
                        plot([timeu(1) timeu(end)],[luy(end) luy(end)],'k:');
                        plot([timeu(1) timeu(end)],[luz(end) luz(end)],'k:'); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                    subplot(4,5,10)  
                           plot(time_l,max(lby,lbz)./max(luy,luz), ...
                                'LineStyle', '-',...
                                'color',color{1},...
                                'LineWidth',1.5,...
                                'MarkerEdgeColor', color{1}, ...
                                'MarkerFaceColor', color{1}, ...
                                'MarkerSize',1.5);
                                title('L_B/L_u'); 
                                xlim([0 t(end)]) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   subplot(4,5,11:12)
                        p01= plot(timeb2,b01,'DisplayName','k=0.1');
                        title('E_m(k^* \leq k)/E_m(k^* > k)-1')
                        hold on;
                        p02= plot(timeb2,b02-1,'DisplayName','k=0.2','LineWidth',1.5);
                        p03= plot(timeb2,b03-1,'DisplayName','k=0.3','LineWidth',1.5);
                        p04= plot(timeb2,b04-1,'DisplayName','k=0.4','LineWidth',1.5);
                        p05= plot(timeb2,b05-1,'DisplayName','k=0.5','LineWidth',1.5);
                        p06= plot(timeb2,b06-1,'DisplayName','k=0.6','LineWidth',1.5);
                        p07= plot(timeb2,b07-1,'DisplayName','k=0.7','LineWidth',1.5);
                        p08= plot(timeb2,b08-1,'DisplayName','k=0.8','LineWidth',1.5);
                        p09= plot(timeb2,b09-1,'DisplayName','k=0.9','LineWidth',1.5);
                        p10= plot(timeb2,b10-1,'DisplayName','k=1.0','LineWidth',1.5);
                        %legend([p01 p02 p03 p04 p05 p06 p07 p08 p09 p10],'Location','west')
                        plot([t(1) t(end)],[0 0],'k:'); 
                        xlim([0 t(end)]) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    subplot(4,5,13:15)
                        p01= plot(timeu2,u01,'DisplayName','k=0.1');
                        title('E_k(k^* \leq k)/E_k(k^* > k)-1')
                        hold on;
                        p02= plot(timeu2,u02-1,'DisplayName','k=0.2','LineWidth',1.5);
                        p03= plot(timeu2,u03-1,'DisplayName','k=0.3','LineWidth',1.5);
                        p04= plot(timeu2,u04-1,'DisplayName','k=0.4','LineWidth',1.5);
                        p05= plot(timeu2,u05-1,'DisplayName','k=0.5','LineWidth',1.5);
                        p06= plot(timeu2,u06-1,'DisplayName','k=0.6','LineWidth',1.5);
                        p07= plot(timeu2,u07-1,'DisplayName','k=0.7','LineWidth',1.5);
                        p08= plot(timeu2,u08-1,'DisplayName','k=0.8','LineWidth',1.5);
                        p09= plot(timeu2,u09-1,'DisplayName','k=0.9','LineWidth',1.5);
                        p10= plot(timeu2,u10-1,'DisplayName','k=1.0','LineWidth',1.5);
                        legend([p01 p02 p03 p04 p05 p06 p07 p08 p09 p10],'Location','eastoutside')
                        plot([t(1) t(end)],[0 0],'k:'); 
                        xlim([0 t(end)]) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

                    subplot(4,5,16)
                        p3= plot(t,w2, ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5);
                            title('Enstrophy');  
                            xlim([0 t(end)])
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
                     subplot(4,5,17)
                        p4= plot(t,hv, ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5);  
                            title('Kinetic Helicity');
                            xlim([0 t(end)])
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                    subplot(4,5,18)
                        p5= plot(t,hc, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5);
                            title('Cross Helicity');  
                            xlim([0 t(end)])
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
                     subplot(4,5,19)
                        p6= plot(t,hm, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5);
                            title('Magnetic Helicity');  
                            xlim([0 t(end)])
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    subplot(4,5,20)
                        p6= plot(t,hm - hv, ...
                        'LineStyle', '-',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',1.5);
                        title('Residual Helicity');  
                        xlim([0 t(end)])
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            color{1} = [0.83,0,0.17]; color{2} = [1,0.67,0.33]; color{3} = [0.33,0.33,1.0];   
            color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5];
            color{6}=[0 252 17]/255;
            color{7}=[16 119 37]/255;
            color{8}=[40 182 163]/255;
            color{9}=[40 59 182]/255;
            color{10}=[130 40 182]/255;

            if RE1_PLOTS ==1 
            if Ree(j) == 1 && Rmm(j) >= 1
                count_ree_1 = count_ree_1 + 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hFig = figure(3000);
                set(hFig, 'Position', [100, 50, 1100, 800]);  
                subplot(4,2,count_to_4 + (ii-1)*2)
                    pre1(count_ree_1)= loglog(t,em, ...
                            'LineStyle', '-',...
                            'color',color{count_ree_1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{count_ree_1}, ...
                            'MarkerFaceColor', color{count_ree_1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName',[ 'R_m=', num2str(Rmm(j))]); 
                if ii == 2 && jj ==1 
                    ylabel('u_{III}', 'FontWeight', 'Bold')
                elseif ii ==1 && jj ==1 
                    ylabel('u_{abc}', 'FontWeight', 'Bold')
                    title(['[',num2str(box(jj)),',2,1]']);
                elseif count_to_4 + (ii-1)*2 == 1 || count_to_4 + (ii-1)*2 == 2
                    title(['[',num2str(box(jj)),',2,1]']);      
                end
%                 if jj == 2
%                     xlabel('Time')
%                 end
                ylim([1e-2 1e2])
                        hold on;
                    pre2(count_ree_1)= loglog(t,ev, ...
                            'LineStyle', '-.',...
                            'color',color{count_ree_1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{count_ree_1}, ...
                            'MarkerFaceColor', color{count_ree_1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','E_k');
                        if ii==2 && jj==2
                        legend([pre1(1:end)],'Location','southeast');%, 'Orientation','Horizontal')
                        legend('boxoff')
                        end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                hFig = figure(3000);
                set(hFig, 'Position', [100, 50, 1100, 800]);  
                subplot(4,2,2 + count_to_4 + (ii-1)*2)
                    semilogx(time_l,max(lby,lbz)/6.28, ...
                            'LineStyle', '-',...
                            'color',color{count_ree_1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{count_ree_1}, ...
                            'MarkerFaceColor', color{count_ree_1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^B_y');  
                        hold on;
                    semilogx(time_l,max(luy,luz)/6.28, ...
                            'LineStyle', '--',...
                            'color',color{count_ree_1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{count_ree_1}, ...
                            'MarkerFaceColor', color{count_ree_1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','L^B_y');      
                if ii == 2 && jj ==1 
                    ylabel('u_{III}', 'FontWeight', 'Bold')
                elseif ii ==1 && jj ==1 
                    ylabel('u_{abc}', 'FontWeight', 'Bold')
%                     title(['[',num2str(box(jj)),',2,1]']);
%                 elseif ii ==1 && jj ==2 
%                     title(['[',num2str(box(jj)),',2,1]']);      
                end
                if  count_to_4 + (ii-1)*2 == 5 || count_to_4 + (ii-1)*2 == 6
                    xlabel('Time')
                end
                ylim([0.1 box(jj)])
                if ii==2 && jj==2
                    l=legend([pre1(1:end)],'Location','southeast', 'Orientation','Horizontal');
                    set(l, 'Position',[0.4, 0.01, .25, .03])
                end
            end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                %%%%%%%%%%%%%%%%%
                % PLOT Lu vs LB %
                %%%%%%%%%%%%%%%%%      
            
            rows = 4;
            if mod(size(Ree,1),rows)==0
                cols = fix(size(Ree,1)/rows);
            else 
                cols = fix(size(Ree,1)/rows) + 1;
            end
%             if run(j) < 100
                colorr{1}=[0 0 0];
                colorr{2}=[153,0,0]/255;
                colorr{3}=[0,102,0]/255;
%             else
%                 colorr{1}=[0 0 0];
%                 colorr{2}=[255,102,102]/255;
%                 colorr{3}=[102,255,102]/255;
%             end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % Plot all available 
                if LENGTH_3D_PLOTS ==1
                hFig = figure(2000 + box(jj) + ii);
                set(hFig, 'Position', [100, 50, 1100, 800]); 
                 
                subplot(cols,rows, j)
                    plot3(lby/(6.28),luy/(6.28),time_l, 'o',...
                        'LineStyle', 'none',...
                        'color',colorr{2},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', colorr{2}, ...
                        'MarkerFaceColor', colorr{2}, ...
                        'MarkerSize',1.5 );
                    grid on;
                    grid minor;
                    hold on;
                    plot3(lbz/(6.28),luz/(6.28),time_l, 'o',...
                        'LineStyle', 'none',...
                        'color',colorr{3},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', colorr{3}, ...
                        'MarkerFaceColor', colorr{3}, ...
                        'MarkerSize',1.5 );
                    plot3(lb2/(6.28),lu2/(6.28),time_l2, 'o',...
                        'LineStyle', 'none',...
                        'color',colorr{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', colorr{1}, ...
                        'MarkerFaceColor', colorr{1}, ...
                        'MarkerSize',1.5 );
                    
                    xlim([0 box(jj)])
                    ylim([0 2])
                    title(['R_e=', num2str(Ree(j)),' R_m=', num2str(Rmm(j))]);
                    if j> size(Ree,1) - rows
                        ylabel(['L_u/2\pi'], 'FontSize', 12);   
                        xlabel(['L_B/2\pi'], 'FontSize', 12);
                    end
                    if mod(j,rows) == 1    
                        zlabel('Time');
                    end
                    plot3(lbz/(6.28),luz/(6.28),zeros(size(timeb)),'k:');
                 end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % Plot Re=1&& Rmm(j) >= 1
                if LENGTH_3D_RE1_PLOTS ==1
                hFig = figure(2050 + box(jj) + ii);
                set(hFig, 'Position', [100, 50, 1200, 200]); 
                rows = 5;
                cols = 1;
%                 if mod(size(Ree(Ree==1 && Rmm >= 1),1),rows)==0
%                     cols = fix(size(Ree(Ree==1),1)/rows);
%                 else 
%                     cols = fix(size(Ree(Ree==1),1)/rows) + 1;
%                 end
%                 
                if Ree(j)==1&& Rmm(j) >= 1
                    count_ree1 = count_ree1 + 1;
                    subplot(cols,rows, count_ree1)
                        plot3(lby/(6.28),luy/(6.28),time_l, 'o',...
                            'LineStyle', 'none',...
                            'color',colorr{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', colorr{2}, ...
                            'MarkerFaceColor', colorr{2}, ...
                            'MarkerSize',1.5 );
                        grid on;
                        grid minor;
                        hold on;
                        plot3(lbz/(6.28),luz/(6.28),time_l, 'o',...
                            'LineStyle', 'none',...
                            'color',colorr{3},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', colorr{3}, ...
                            'MarkerFaceColor', colorr{3}, ...
                            'MarkerSize',1.5 );
                        plot3(lb2/(6.28),lu2/(6.28),time_l2, 'o',...
                            'LineStyle', 'none',...
                            'color',colorr{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', colorr{1}, ...
                            'MarkerFaceColor', colorr{1}, ...
                            'MarkerSize',1.5 );

                        xlim([0 box(jj)])
                        ylim([0 2])
                        title(['R_e=', num2str(Ree(j)),' R_m=', num2str(Rmm(j))]);
                    
                        ylabel(['L_u/2\pi'], 'FontSize', 12);   
                        xlabel(['L_B/2\pi'], 'FontSize', 12);
     
                    if count_ree1  == 1    
                        zlabel('Time');
                    end
                    plot3(lbz/(6.28),luz/(6.28),zeros(size(timeb)),'k:');   
                end
                end
                %%%%%%%%%%%%%%
                % SAVE EM/EK %
                %%%%%%%%%%%%%%
                emf(j) = em(end);
                evf(j) = ev(end);     
                lbf(j) = max(lby(end),lbz(end));
                lvf(j) = max(luy(end),luz(end));

                if run(j) < 100
                    count1 = count1 + 1;
                    ratio1(count1) = em(end)./ev(end); 
                    ratio_l1(count1) = lbf(j)./lvf(j);
                    Ree1(count1) = Ree(j);
                    Rmm1(count1) = Rmm(j);
                else
                    count2 = count2 + 1;
                    ratio2(count2) = em(end)./ev(end); 
                    ratio_l2(count2) = lbf(j)./lvf(j);
                    Ree2(count2) = Ree(j);
                    Rmm2(count2) = Rmm(j);
                end
        end
        end  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % In order to compare the theoretical results with the sims
        % we need to scale them accordingly: 
        
        if ENERGY_BOUND_PLOTS ==1
            color{5} = [255 179 9]/255;     color{2} = [173 95 252]/255; 
            color{3} = [0.33,0.33,1.0];     color{4} = [0,0.5,0]; 
            color{1} = [39 235 242]/255;
            
            hFig = figure(5000);
            set(hFig, 'Position', [100, 50, 1200, 600]); 
            re_here = 0.1:0.9:100;
            rm_here = [1 10 20 50 100];
            evf_here = 0.05:0.05:2;
            nu = 1./re_here;                               
            eta= 1./rm_here;     
            G = real(sqrt(1./(2*evf_here)));
            [NU,GG] = meshgrid(nu,G);
            FV = NU;
            if ii == 1
                FV = sqrt(3);      % F times nu
            else 
                FV = NU.* sqrt(1-6.*NU.^2)/(2.*NU.^2);
            end
            
            if jj==1  % to avoid repetition
                if ii==1
                    lower = real(sqrt(NU.*((FV.*GG).^2 - FV.*GG) -1));
                else
                    lower = real(sqrt(0.4.*NU.*(FV.*GG).^2 )- 2/sqrt(8).*(NU).^2./(sqrt(1-6*NU.^2)).*FV.*GG -1);
                end
                
                for ind = 1:size(rm_here,2)
                    pm_here = rm_here(ind)./re_here;
                    [PM, GG] = meshgrid(pm_here, G);
                    UP = FV.*GG;
                    %upper = PM.*(GG -1);
                    for ind1=1:size(GG,1)
                        for ind2=1:size(GG,2)
                                upper(ind1,ind2) = PM(ind1,ind2)*(UP(ind1,ind2)-1); 
                        end
                    end
               
                subplot(1,2,ii)
                    s3(ind) =surf(1./(2.*(FV.*GG).^2),1./NU,upper, 'DisplayName',['R_m=', num2str(rm_here(ind))],'FaceColor',color{ind},'EdgeColor','black');
                    set(s3(ind),'FaceAlpha',0.6- ind*0.05,'EdgeAlpha',0.5- ind*0.05);
                    hold on;
                    view(60,10)   
                    xlabel('E_k')
                    ylabel('R_e') 
                    zlabel('|E_m/E_k|')
                    %shading interp
                    if ii == 1
                        title(['u_{ABC}'])
                    else
                        title(['u_{III}'])
                    end
                    zlim([0 10])
                    set(get(s3(ind),'Parent'),'YScale','log')%,'ZScale','log');
                    
                end
                
                s1(ii) = surf(1./(2.*(FV.*GG).^2),1./NU,lower,'Marker','o','MarkerSize',1.0,'FaceColor','red','EdgeColor','black', 'DisplayName','Lower Limit');
                if ii ==1 
                    l = legend([s1,s3]); 
                    set(l,'Location','best')
                    set(l,'Position',[0.42 0.82 0.1 0.02])
                end
             end
            clear NU GG lower upper FV UP Fv
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Overplot the results from sims
         
            nu_res = 1./transpose(Ree);
            eta_res= 1./transpose(Rmm);
%             G_res = real(sqrt(1./(2*evf)));
            sz = 30;
            color{5} = [255 179 9]/255;     color{2} = [173 95 252]/255; 
            color{3} = [0.33,0.33,1.0];     color{4} = [0,0.5,0]; 
            color{1} = [39 235 242]/255;
            for ind1=1:size(nu_res,2)
                if Rmm(ind1) < rm_here(1)
                    c(ind1,1) = 39/255.0;
                    c(ind1,2) = 235/255.0;
                    c(ind1,3) = 235/255.0;
                elseif Rmm(ind1) < rm_here(2)
                    c(ind1,1) = 173/255.0;
                    c(ind1,2) = 95/255.0;
                    c(ind1,3) = 252/255.0;
                    
                elseif Rmm(ind1) < rm_here(3)
                    c(ind1,1) = 0.33;
                    c(ind1,2) = 0.33;
                    c(ind1,3) = 1.0;

                elseif Rmm(ind1) < rm_here(4)
                    c(ind1,1) = 0.0;
                    c(ind1,2) = 0.5;
                    c(ind1,3) = 0.0;                    
                else
                    c(ind1,1) = 255/255.0;
                    c(ind1,2) = 179/255.0;
                    c(ind1,3) = 9/255.0;
  
                end
            end
            %scatter(x,y,sz,c,'filled')
            subplot(1,2,ii) 
                if ii == 1
                    scatter3(evf,1./nu_res,emf./evf ,sz,c,'filled','MarkerEdgeColor','k')
                end
                if ii == 2
                    scatter3(evf,1./nu_res,emf./evf ,sz,c,'filled','MarkerEdgeColor','k') 
                end 
            clear c
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if LENGTH_3D_PLOTS ==1
        figure(2000 + box(jj) + ii);
        if ii ==1 
            suptitle(['u_{abc} in [',num2str(box(jj)),',2,1]']);
        else
            suptitle(['u_{III} in [',num2str(box(jj)),',2,1]']);
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if LENGTH_3D_RE1_PLOTS ==1
        figure(2050 + box(jj) + ii);
        if ii ==1 
            suptitle(['u_{abc} in [',num2str(box(jj)),',2,1]']);
        else
            suptitle(['u_{III} in [',num2str(box(jj)),',2,1]']);
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if RATIO_PLOT ==1
        hFig = figure(1000 + jj);
            set(hFig, 'Position', [100, 50, 500, 500]); 
            if ii == 1
                color{1}=[0.83 0.17 0.17];
            else
                color{1}=[0 0 0];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(2,3,1)
            if exist('Rmm1','var') == 1
                p11 = semilogx(Rmm1, ratio1, 'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );%,...
                        %'DisplayName',['[',num2str(box(jj)),',2,1]']
                        hold on;
            end
            if exist('Rmm2','var') == 1
                p11 = semilogx(Rmm2, ratio2, '^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );
                    hold on;
                        %xlim([0.2 10])
                        %ylim([0 6])
            end
                        ylabel(['E_m/E_k']);  
                        xlabel('R_m')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            subplot(2,3,2)
            if exist('Ree1','var') == 1
                p11 = semilogx(Ree1, ratio1, 'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );%,...
                    hold on;
            end
             if exist('Ree2','var') == 1           
               p11 = semilogx(Ree2, ratio2, '^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );
                    hold on;
            end
                        title(['[',num2str(box(jj)),',2,1]']);
                        xlabel('R_e')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(2,3,3)
            if exist('Ree1','var') == 1
                p13 = semilogx(Rmm1./Ree1,ratio1,'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5);
                          hold on;
            end
             if exist('Ree2','var') == 1                  
                p13 = semilogx(Rmm2./Ree2,ratio2,'^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5); 
            end
                        xlabel('Pr_m')
                        %title('E_m')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             subplot(2,3,4)
             if exist('Ree1','var') == 1
                p11 = semilogx(Rmm1, ratio_l1, 'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );
                        hold on;
             end
             if exist('Ree2','var') == 1                    
                p11 = semilogx(Rmm2, ratio_l2, '^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );  
                        hold on;
             end
                        ylabel(['L_B/L_u']);  
                        xlabel('R_m')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(2,3,5)
            if exist('Ree1','var') == 1
                p11 = semilogx(Ree1, ratio_l1, 'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );%,...
                        hold on;
            end
            if exist('Ree2','var') == 1           
                p11 = semilogx(Ree2, ratio_l2, '^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5 );
                           hold on;
            end
                        xlabel('R_e')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(2,3,6)
            if exist('Ree1','var') == 1
                p13 = semilogx(Rmm1./Ree1,ratio_l1,'o',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5);

                        hold on;       
            end
            if exist('Ree2','var') == 1                    
                p13 = semilogx(Rmm2./Ree2,ratio_l2,'^',...
                        'LineStyle', 'none',...
                        'color',color{1},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{1}, ...
                        'MarkerFaceColor', color{1}, ...
                        'MarkerSize',4.5);
                        hold on;
            end
                         xlabel('Pr_m')

            end
            %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%           
            hFig = figure(1000);
            set(hFig, 'Position', [100, 50, 900, 600]); 
            color{1}=[0.83 0.17 0.17];
            color{2}=[0 0 0];

            subplot(1,2,jj)

             if exist('Ree1','var') == 1
                p = loglog(Rmm1./Ree1, ratio1./ratio_l1, 'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5 );%,...
                        %'DisplayName',['[',num2str(box(jj)),',2,1]']

                        hold on;
             end
            if exist('Ree2','var') == 1                   
                p = loglog(Rmm2./Ree2, ratio2./ratio_l2, '^',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5 );
                     hold on;
             end
                        xlabel('Pr_m')
                        if jj == 1
                            ylabel(['E_m/E_k /  L_B/L_u']);  
                        end
                        title(['[',num2str(box(jj)),',2,1]']);
                        if ii == 2
                            legend('u_{abc} saturated','u_{abc} in progress', 'u_{III} saturated','u_{III} in progress','Location','southeast');
                            legend('boxoff');
                        end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% KINEMATIC R %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if KINEMATIC_CF_PLOTS == 1
            rep3=[rep2{ii},'/results/box_', num2str(box(jj))];
            file_infokr=[rep3, '/info'];
            info_filekr=importdata(file_infokr);
            infokr=info_filekr;
            tblAkr = table(infokr(:,1),infokr(:,2), infokr(:,3));
            % Sort the rows of the table based on Rm
            tblBkr = sortrows(tblAkr,3); 
            runkr = tblBkr{1:end,1}; 
            caskr = tblBkr{1:end,2}; 
            Rmmkr = tblBkr{1:end,3};
            for jjj=1:size(Rmmkr,1)
                file_timekr=[rep2{ii},'/kinematicOutput_box_', num2str(box(jj)), '_' ,num2str(runkr(jjj)), '/result',num2str(caskr(jjj)), '/timevar2'];
                full_filekr=importdata(file_timekr);
                timevarkr=full_filekr;%.data;
                lbxkr(jjj)  =timevarkr(end,5);
                lbykr(jjj)  =timevarkr(end,6);
                lbzkr(jjj)  =timevarkr(end,7);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hFig = figure(4000);
                color{1}=[0.83 0.17 0.17];
                color{2}=[0 0 0];
                set(hFig, 'Position', [100, 50, 1100, 800]);  
                subplot(2,2,count_to_4)
                    p_nl = semilogx(Rmm,lbf./6.28, 'o', ...
                            'LineStyle', 'none',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',3.5,...
                            'DisplayName','nonlinear'); 
                        hold on;
                    p_kr = semilogx(Rmmkr,max(lbykr,lbzkr)/6.28, 'o', ...
                            'LineStyle', 'none',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',3.5,...
                            'DisplayName','kinematic');   
                if ii == 2 && jj ==1 
                    ylabel('u_{III}', 'FontWeight', 'Bold')
                elseif ii ==1 && jj ==1 
                    ylabel('u_{abc}', 'FontWeight', 'Bold')
                    title(['[',num2str(box(jj)),',2,1]']);
                elseif ii ==1 && jj ==2 
                    title(['[',num2str(box(jj)),',2,1]']);      
                end
                if ii == 2
                    xlabel('R_m')
                end
                if ii ==2 && jj ==2 
                    legend([p_nl, p_kr])
                end
                ylim([0,box(jj) ])
            end    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    end
           
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hFig = figure(1000);
  Prm = sort(Rmm2./Ree2);
 subplot(1,2,1)
 plot(Prm, 0.3*Prm.^0.33, ':k')
 ylim([1e-4 10])
 subplot(1,2,2)
 plot(Prm, 0.3*Prm.^0.5, ':k')
  ylim([1e-4 10])
        
      