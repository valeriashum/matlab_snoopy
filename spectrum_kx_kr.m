close all;
clear all;
rep1=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc/', '/store/ASTRO/vs391/kinematic_dynamo/u_iii/'}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET ?_PLOTS=1 to get all plots     %%
% SET ?_PLOTS=0 to hide state plots  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%
WATERFALL_PLOTS = 0;    %
ALL_SPECTRUM_PLOTS = 0; %
SPECTRUM_32_PLOTS = 0;  %
MAX_KX_PLOTS = 1;       %
%%%%%%%%%%%%%%%%%%%%%%%%%
box=[8 10 16 20 32];

rmc_iii=[0.627 0.54  0.4 0.35 0.27 ]; 
rmc_abc=[4.275000e-01 0.3385 0.2695  0.2354 0.1753]; %abc

rows = 1;
if mod(size(box,1),rows)==0
    cols = fix(size(box,1)/rows);
else 
    cols = fix(size(box,1)/rows) + 1;
end

color{1}=[0.83 0.17 0.17];  %abc
color{2}=[0 0 0];           %iii

color{1}=[255 0 0]/255;
color{2}=[216 48 48]/255;
color{3}=[255 108 9]/255;
color{4}=[250 213 24]/255;
color{5}=[210 252 0]/255;
color{6}=[0 252 17]/255;
color{7}=[16 119 37]/255;
color{8}=[40 182 163]/255;
color{9}=[40 59 182]/255;
color{10}=[130 40 182]/255;
color{11}=[255 0 0]/255;
color{12}=[216 48 48]/255;
color{13}=[255 108 9]/255;
color{14}=[250 213 24]/255;
color{15}=[210 252 0]/255;
color{16}=[0 252 17]/255;
color{17}=[16 119 37]/255;
color{18}=[40 182 163]/255;
color{19}=[40 59 182]/255;
color{20}=[130 40 182]/255;
color{21}=[255 0 0]/255;
color{22}=[216 48 48]/255;
color{23}=[255 108 9]/255;
color{24}=[250 213 24]/255;
color{25}=[210 252 0]/255;
color{26}=[0 252 17]/255;
color{27}=[16 119 37]/255;
color{28}=[40 182 163]/255;
color{29}=[40 59 182]/255;
color{30}=[130 40 182]/255;
color{31}=[255 0 0]/255;
color{32}=[216 48 48]/255;
color{33}=[255 108 9]/255;
color{34}=[250 213 24]/255;
color{35}=[210 252 0]/255;
color{36}=[0 252 17]/255;
color{37}=[16 119 37]/255;
color{38}=[40 182 163]/255;
color{39}=[40 59 182]/255;

for ii=1:size(rep1,2)
    clear tblA tblB run cas Rmm max_kx max_rm
    count = 0;
    
    for jj=1:size(box,2)
        clear ksq emx emy emz p1
        count = count + 1;
        rep3=[rep1{ii},'/results/box_', num2str(box(jj))];
        file_info=[rep3, '/info'];
        full_file=importdata(file_info);
        info=full_file;
        tblA = table(info(:,1),info(:,2), info(:,3));
        % Sort the rows of the table based on Re, then Rm
        tblB = sortrows(tblA,[3]); 
        run = tblB{1:end,1}; 
        cas = tblB{1:end,2}; 
        Rmm  = tblB{1:end,3};
        
        if ii == 2 
            m = 1.0/box(jj);
            Rmm = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Rmm;
        end     
        
        count_rmm = 0;
        for j=1:size(Rmm,1)   
            if ii==1 && Rmm(j) < rmc_abc(jj)
                fprintf('R_m out of range: R_m = %f < R_m^c=%f,\n',Rmm(j),rmc_abc(jj))
            elseif ii==2 && Rmm(j) < rmc_iii(jj)
                fprintf('R_m out of range\n')
            else
                clear x y z
                count_rmm = count_rmm + 1;
                file_time=[rep1{ii},'/kinematicOutput_box_', num2str(box(jj)), '_' ,num2str(run(j)), '/result',num2str(cas(j)), '/spectrum_kx.dat'];
                full_file=importdata(file_time);
                timevar=full_file;
                tblA = table(timevar(:,1),timevar(:,2), timevar(:,3), timevar(:,4));
                ksq = tblA{1:end,1}/box(jj)^2;
                emx = tblA{1:end,2};
                emy = tblA{1:end,3};
                emz = tblA{1:end,4};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if MAX_KX_PLOTS == 1
                [M,I] = max(emx + emy + emz);
                max_kx(count_rmm) = ksq(I);
                max_rm(count_rmm) = Rmm(j);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if WATERFALL_PLOTS == 1
                hFig = figure(ii);
                set(hFig, 'Position', [100, 50, 1200, 900]);  

                    z(:,1) = (emx + emy + emz)/max(emx + emy + emz);
                    x = [ Rmm(j)];
                    y = sqrt(ksq);
                subplot(3,2, jj)
                w = waterfall(x,y,z);
                    %zlim([0.95 1])
                    w.Marker = 'o'; 
                    w.MarkerSize = 3.5;
                    w.MarkerFaceColor = color{j};
                    w.MarkerEdgeColor = color{j};
                    w.LineStyle = ':';
                    w.FaceColor = color{j};
                    hold on;
                    grid on;
                    grid minor;
                    set(gca, 'YScale','log')
                    title(['[',num2str(box(jj)),',2,1]']);
                    xlabel('R_m')
                    ylabel('|k_x|')
                    zlabel('Normalised E_m')
                    view(60, 10);
                end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ALL_SPECTRUM_PLOTS == 1
                hFig = figure(10 + ii);
                set(hFig, 'Position', [100, 50, 1200, 900]);  
                subplot(1,5, jj)
                loglog(sqrt(ksq),(emx + emy + emz)/max(emx + emy + emz),'o',...
                            'LineStyle', 'none',... 
                            'color',color{j},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{j}, ...
                            'MarkerFaceColor', color{j}, ...
                            'MarkerSize',3.5,...
                            'DisplayName',['R_m=',num2str(Rmm(j))]); 
               title(['[',num2str(box(jj)),',2,1]']);
                    xlabel('|k_x|')
                    ylabel('Normalised E_m')  
                    hold on
                end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if SPECTRUM_32_PLOTS ==1 
                hFig = figure(100);
                set(hFig, 'Position', [100, 50, 1200, 900]);  
               
                if jj == 5 
                    ax(ii) = subplot(3,2, [ii,ii+2])
                    p1(count_rmm) = loglog(sqrt(ksq),(emx + emy + emz)/max(emx + emy + emz),'o',...
                            'LineStyle', 'none',... 
                            'color',color{j},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{j}, ...
                            'MarkerFaceColor', color{j}, ...
                            'MarkerSize',3.5,...
                            'DisplayName',['R_m=',num2str(Rmm(j))]); 
                        if ii==1
                            title(['u_{abc}']);
                            ylabel('Normalised E_m')
                        else 
                            title(['u_{III}']);
                        end

                        hold on 
                        legendInfo{1,count_rmm} = {['R_m=',num2str(Rmm(j))]};
                    ax(2 + ii) = subplot(3,2, 4 + ii)
                    loglog(sqrt(ksq),(emx + emy + emz)/max(emx + emy + emz),'o',...
                            'LineStyle', ':',... 
                            'color',color{j},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{j}, ...
                            'MarkerFaceColor', color{j}, ...
                            'MarkerSize',3.5,...
                            'DisplayName',['R_m=',num2str(Rmm(j))]);  
                        hold on
                        xlabel('|k_x|')
                        ylim([0.95 1.0])
                        xlim([0.01 100])
                    linkaxes([ax(1:end)],'x')
                end
               end   
            end
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if SPECTRUM_32_PLOTS ==1 
            if jj == 5
                 figure(100)
                 l = legend([p1(1:end)],'Location','southwest', 'Orientation','vertical');
                 plot([1 1],[1e-4 1.0], 'k--');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if MAX_KX_PLOTS == 1 
        hFig = figure(20);
        set(hFig, 'Position', [100, 50, 600, 500]);  
        subplot(1,2,ii)
        p2(jj) = plot(max_rm, sqrt(max_kx),'o',...
                            'LineStyle', 'none',... 
                            'color',color{jj*2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{jj*2}, ...
                            'MarkerFaceColor', color{jj*2}, ...
                            'MarkerSize',6.5,...
                            'DisplayName',['[',num2str(box(jj)),',2,1]']); 
                        hold on;
                        xlabel('R_m')
                        if ii == 1
                            ylabel('Dominant k_x')
                        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MAX_KX_PLOTS == 1 
        figure(20)
        if ii == 2
            l = legend([p2(1:end)]);%,'Location','southwest', 'Orientation','vertical');
        end
    end
end