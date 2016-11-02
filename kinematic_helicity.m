clear all;
close all;
rep1=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc', '/store/ASTRO/vs391/kinematic_dynamo/u_iii'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET ALL_PLOTS=1 to get all plots     %%
% SET ALL_PLOTS=0 to hide state plots  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
ALL_PLOTS = 1;
box=[8 10 16 20 32];
color{1} = [0.83,0,0.17];       color{2} = [1,0.67,0.33]; 
color{3} = [0.33,0.33,1.0];     color{4} = [0,0.5,0]; 
color{5} = [0.5,0,0.5];

rmc_iii=[0.627 0.54  0.4 0.35 0.27 ]; 
rmc_abc=[4.275000e-01 0.3385 0.2695  0.2354 0.1753];

for ii=1:size(rep1,2)
    clear tblA tblB run cas Rmm mh ch p3
    for jj=1:size(box,2)
        clear mh ch
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
        if ii == 1 
            m = 1.0/box(jj);
            Rmm = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Rmm;
        end
        
        for j=1:size(Rmm,1)    
            file_time=[rep1{ii},'/kinematicOutput_box_', num2str(box(jj)), '_' ,num2str(run(j)), '/result',num2str(cas(j)), '/timevar6'];
            full_file=importdata(file_time);
            timevar=full_file;
            mh(j) = timevar(end,2);
            ch(j) = timevar(end,1);
        end  
        
        hFig = figure(1);
        set(hFig, 'Position', [100, 50, 1100, 900]); 
        subplot(2,2,ii)
            p3(jj)= loglog(Rmm,abs(mh),'o', ...
                        'LineStyle', 'none',...
                        'color',color{jj},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{jj}, ...
                        'MarkerFaceColor', color{jj}, ...
                        'MarkerSize',4.5,...
                        'DisplayName',['[', num2str(box(jj)), ',2,1']); 
                    
                if ii == 1   
                    ylabel('Magnetic Helicity') 
                    title('u_{abc}')
                else
                    title('u_{III}')
                end
                hold on;
                if ii == 1 
                    plot([rmc_abc(jj) rmc_abc(jj)],[1e-15 1e10],':','Color', color{jj});
                else 
                    plot([rmc_iii(jj) rmc_iii(jj)],[1e-15 1e10],':','Color', color{jj});
                end
        subplot(2,2,ii+2)
            p3(jj)= loglog(Rmm,abs(ch),'o', ...
                        'LineStyle', 'none',...
                        'color',color{jj},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{jj}, ...
                        'MarkerFaceColor', color{jj}, ...
                        'MarkerSize',4.5,...
                        'DisplayName',['[', num2str(box(jj)), ',2,1]']); 
                    ylim([1e-25 1])
                if ii == 1    
                    ylabel('Cross Helicity') 
                end
                    xlabel('R_m')
                
                hold on;
                if ii == 1 
                    plot([rmc_abc(jj) rmc_abc(jj)],[1e-25 1],':','Color', color{jj});
                else 
                    plot([rmc_iii(jj) rmc_iii(jj)],[1e-25 1],':','Color', color{jj});
                end
    end
    if ii == 1
         l = legend([p3(1:size(box,2))],'Location','north', 'Orientation','Horizontal')
         set(l, 'Position',[0.40, 0.49, .25, .05])
    end
    
end