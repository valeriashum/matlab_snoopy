% Timevar from the hydro sims
close all
clear all
rep1=[{'/data/novadisk/vs391/hydro/u_abc', '/data/novadisk/vs391/hydro/u_iii'}];
box=[10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET ?_PLOTS=1 to get all plots     %%
% SET ?_PLOTS=0 to hide state plots  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%
WATERFALL_PLOTS = 1;    %
ENERGY_PLOTS = 0;       %
%%%%%%%%%%%%%%%%%%%%%%%%%

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

for ii=1:size(rep1,2)
    clear tblA tblB run cas Ree
    count = 0;
    
    for jj=1:size(box,2)
        clear run cas Ree t ev hv w2 p1
        count = count + 1;
        % Read from info file what Re's we have
        % and the corresponding location
        rep3=[rep1{ii},'/box_', num2str(box(jj))];
        file_info=[rep3, '/info'];
        full_file1=importdata(file_info);
        timevar1=full_file1;%.data;
        tblA = table(timevar1(:,1),timevar1(:,2), timevar1(:,3));
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,3); 
        run = tblB{1:end,1}; 
        cas = tblB{1:end,2}; 
        Ree  = tblB{1:end,3}; 
        
        if ii == 2 
            m = 1.0/box(jj);
            Ree = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Ree;
        end     
      
        count_ree = 0;
        for kk=1:size(Ree,1)
            clear x y z
            count_ree = count_ree + 1;
            file_time=[rep1{ii},'/box_',num2str(box(jj)), '/run_', num2str(run(kk)),'/timevar',num2str(cas(kk)),'.dat'];
            full_file=importdata(file_time);
            timevar=transpose(full_file.data);
            nvar=size(timevar,1);
            var_name=strread(full_file.textdata{1},'%s',nvar);

            % ASSIGN VARIABLE WITH HEADINGS 
            for j=1:nvar
                assignin('base',var_name{j},timevar(j,:)); 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ENERGY_PLOTS == 1
            hFig = figure(1);
            set(hFig, 'Position', [100, 50, 1200, 900]); 
            subplot(3,2,ii)
                p1(ii) = semilogy(t,ev, ...
                    'LineStyle', '-',...
                    'color',color{kk},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color{kk}, ...
                    'MarkerFaceColor', color{kk}, ...
                    'MarkerSize',1.5,...
                    'DisplayName',['R_e=',num2str(Ree(kk))]);
                if ii == 1
                    ylabel('Kinetic Energy');
                end
                if ii == 1
                    title('u_{abc}')
                else 
                    title('u_{III}')
                end
            hold on;        
            subplot(3,2,ii+2)
                plot(t,w2, ...
                    'LineStyle', '-',...
                    'color',color{kk},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color{kk}, ...
                    'MarkerFaceColor', color{kk}, ...
                    'MarkerSize',1.5);
                if ii == 1
                    ylabel('Enstrophy');  
                end
             hold on;  
             subplot(3,2,ii+4)
                plot(t,hv, ...
                    'LineStyle', '-',...
                    'color',color{kk},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color{kk}, ...
                    'MarkerFaceColor', color{kk}, ...
                    'MarkerSize',1.5);  
                if ii == 1
                    ylabel('Kinetic Helicity');
                end   
                xlabel('Time')
             hold on;  
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if WATERFALL_PLOTS == 1
            hFig2 = figure(2);
            set(hFig2, 'Position', [100, 50, 1200, 900]); 
            color{2}=[0 0 0]; 
            for ll=1:3 % for ev hv w2
                subplot(3,2,ii + (ll-1)*2)
                if ll == 1 
                    z(:,1) = ev;
                elseif ll ==2 
                    z(:,1) = abs(hv);
                else
                    z(:,1) = abs(w2);
                end
                x = [ Ree(kk)];
                y = t;
                    w = waterfall(x,y,z);
                        %zlim([0.95 1])
                        w.Marker = 'o'; 
                        w.MarkerSize = 1.5;
                        w.MarkerFaceColor = color{2};
                        w.MarkerEdgeColor = color{2};
                        w.LineStyle = ':';
                        w.LineWidth = 0.1;
                        w.FaceColor = color{2};
                        hold on;
                        %grid on;
                        %grid minor;
                        set(gca, 'ZScale','log')
                        xlim([0 5])
                    if ii == 1 && ll == 1
                        title('u_{abc}')
                    elseif  ii == 2 && ll == 1
                        title('u_{III}')
                    end
                        xlabel('R_e')
                        ylabel('Time')
                        if ll == 1 && ii == 1
                             zlabel('Kinetic Energy')
                        elseif ll ==2 && ii == 1
                            zlabel('Enstrophy')
                        elseif  ll ==3 && ii == 1
                            zlabel('Kinetic Helicity')
                        end
                       colormap( [ 229 204 255]/255 );
                        view(110, 10);
            end
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if ENERGY_PLOTS == 1
            l = legend([p1(1:end)]);%,'Location','southwest', 'Orientation','vertical');
         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
