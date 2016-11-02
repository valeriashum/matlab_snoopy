close all;
clear all;
rep1=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc/','/store/ASTRO/vs391/kinematic_dynamo/u_iii/'}];
box=[32];
colorr{1}=[0.83 0.17 0.17];  %abc
colorr{2}=[0 0 0];           %iii

color{1}=[255 0 0]/255;
color{2}=[216 48 48]/255;
color{3}=[255 108 9]/255;
color{4}=[250 213 24]/255;
color{5}=[210 252 0]/255;
color{6}=[0 252 17]/255;
color{7}=[16 119 37]/255;
color{8}=[130 40 182]/255;
color{9}=[12 232 247]/255;
color{10}=[40 182 163]/255;
color{11}=[40 59 182]/255;

for ii=1:size(rep1,2)
    clear tblA tblB tblA1 tblB1 run cas Rmm 

    for jj=1:size(box,2)
        clear tblA tblB tblA1 tblB1 run cas Rmm 
        
        rep3=[rep1{ii},'/results/box_', num2str(box(jj))];
        file_info=[rep3, '/info'];
        full_file=importdata(file_info);
        info=full_file;
        tblA1 = table(info(:,1),info(:,2), info(:,3));
        % Sort the rows of the table based on Re, then Rm
        tblB1 = sortrows(tblA1,[3]); 
        run = tblB1{1:end,1}; 
        cas = tblB1{1:end,2}; 
        Rmm = tblB1{1:end,3};
       
        
        rows = 4;
        if mod(size(Rmm,1),rows)==0
            cols = fix(size(Rmm,1)/rows);
        else 
            cols = fix(size(Rmm,1)/rows) + 1;
        end
        count_l = 1;
        count_s = 1;
        for j=1:size(Rmm,1)    
            clear tblA t Ratio timevar
            file_time2=[rep1{ii},'/kinematicOutput_box_', num2str(box(jj)), '_' ,num2str(run(j)), '/result',num2str(cas(j)), '/b_measures_copy_vs_t'];
%             delimiterIn = '+e';
%             full_file2=importdata(file_time2,delimiterIn);
            full_file2=dlmread(file_time2);
            timevar=full_file2;
            tblA = table(timevar(:,1),timevar(:,2),timevar(:,3), ...
                        timevar(:,4),timevar(:,5),timevar(:,6), ...
                        timevar(:,7),timevar(:,8),timevar(:,9), ...
                        timevar(:,10),timevar(:,11),timevar(:,12));
            t  = tblA{1:end,1};
            lb = tblA{1:end,2};
                      
            for ind=1:10
                clear Ratio;
                %color{ind}=rand(1,3);
                Ratio = tblA{:,2+ind};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hFig = figure(ii);
                set(hFig, 'Position', [0, 0, 1200, 900]); 
                set(gca, 'FontSize', 14)
                subplot(cols,rows, j)
                    semilogy(t,Ratio,'o',...
                            'LineStyle', 'none',... 
                            'color',color{ind},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{ind}, ...
                            'MarkerFaceColor', color{ind}, ...
                            'MarkerSize',3.5);     
                     hold on;
                     xlabel('Time')
                     if mod(j,cols) == 1
                        ylabel('E_m^{large}/E_m^{small}')
                     end
                     if j > size(Rmm,1) - cols + 1 
                        xlim([t(1) t(end)])
                     end
                     title(['R_m=',num2str(Rmm(j))])
            end
            if j ==1  
                l = legend({'k_x=0.1','k_x=0.2',...
                   'k_x=0.3','k_x=0.4',...
                   'k_x=0.5','k_x=0.6',...
                   'k_x=0.7','k_x=0.8',...
                   'k_x=0.9','k_x=1.0'},...
                   'Location','eastoutside','Orientation', 'Horizontal');
                set(l, 'Position',[0.10, 0.95, .8, .05])
            end
            plot([t(1) t(end)],[1 1],'k:')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear Ratio z x y
            Ratio = tblA{:,8};
            z(:,1) = Ratio;
            x = [ Rmm(j)];
            y = t;
            

            hFig = figure(100 + ii);
            set(hFig, 'Position', [0, 0, 1200, 900]); 
            set(gca, 'FontSize', 14);
                w = waterfall(x,y,z);
                    hold on;
                    ylim([0 200])
                    w.Marker = 'o'; 
                    w.MarkerSize = 3.5;
                    w.MarkerFaceColor = colorr{ii};
                    w.MarkerEdgeColor = colorr{ii};
                    w.LineStyle = ':';
                    w.FaceColor = colorr{ii};
                    hold on;
                    grid on;
                    grid minor;
                    set(gca, 'ZScale','log')
                    view(100, 10);
                
                ylabel('R_m')
                ylabel('Time')
                zlabel('E_m^{large}/E_m^{small}')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the rate of E^L/E_S growing / falling in time
            % for the scale k=k^u/10
            clear Ratio 
            Ratio = tblA{:,3+1};
           
            if Ratio(end) > 1.0
                rate_l(count_l) = log(Ratio(end)/Ratio(end-5)) / (t(end) - t(end-5));
                rm_l(count_l) = Rmm(j); 
                count_l = count_l + 1;  
            else
                rate_s(count_s) = log(Ratio(end)/Ratio(end-5)) / (t(end) - t(end-5));
                rm_s(count_s) = Rmm(j); 
                count_s = count_s + 1;   
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rep4=[rep1{ii},'/results/box_', num2str(box(jj))];
        file_info4=[rep4, '/result_summary'];
        full_file4=importdata(file_info4);
        timevar4=transpose(full_file4.data);
        nvar=size(timevar4,1);
        var_name=strread(full_file4.textdata{1},'%s',nvar);

        Rm = timevar4(2,:); 
        Gr = timevar4(3,:); 
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure(100)
%         l = legend([p1(1:size(Rmm,2))],'Location','north', 'Orientation','Horizontal')
%         set(l, 'Position',[0.40, 0.49, .25, .05])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig = figure(10 + ii);
        set(hFig, 'Position', [0, 0, 1200, 900]); 
        set(gca, 'FontSize', 14)
        if count_l > 1
            plot(rm_l,rate_l,'o',...
                'LineStyle', 'none',... 
                'color',colorr{1},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', colorr{1}, ...
                'MarkerFaceColor', colorr{1}, ...
                'MarkerSize',3.5);
            hold on;
        end
        if count_s > 1
        plot(rm_s,-rate_s,'o',...
            'LineStyle', 'none',... 
            'color',colorr{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', colorr{2}, ...
            'MarkerFaceColor', colorr{2}, ...
            'MarkerSize',3.5);
            hold on
        end
        plot(Rm,Gr,'^',...
            'LineStyle', 'none',... 
            'color',colorr{2},...
            'LineWidth',1.5,...
            'MarkerEdgeColor', colorr{2}, ...
            'MarkerFaceColor', colorr{2}, ...
            'MarkerSize',3.5);
            hold on
        xlabel('R_m')
        ylabel('Growth rate ratio')
    end
end

