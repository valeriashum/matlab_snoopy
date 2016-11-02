close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep1=[{'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/box_', '/store/ASTRO/vs391/nonlinear_dynamo/u_iii/box_'}];
rep2=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc'     , '/store/ASTRO/vs391/kinematic_dynamo/u_iii'}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box=[10 16 32];
color{1}=[0.83 0.17 0.17];
color{2}=[0 0 0];
color{3}=[160 160 160]/255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;
for ii=1:size(rep1,2)
    clear tblA tblB run cas 
    for jj=1:size(box,2)
        clear em ev yneg ypos  Ek nuinv etainv  
        rep=[rep1{ii}, num2str(box(jj))];
        file_info=[rep, '/info'];
        full_file=importdata(file_info);
        info=full_file;
        tblA = table(info(:,1),info(:,2), info(:,3), info(:,4));
        % Sort the rows of the table based on Re, then Rm
        tblB = sortrows(tblA,[4,3]); 
        run = tblB{1:end,1}; 
        cas = tblB{1:end,2}; 
        etainv = tblB{1:end,3};  
        nuinv  = tblB{1:end,4}; 
        
        for j=1:size(nuinv,1)
            clear em ev t upper lower Pm k2
            if cas(j) ~= 0 && run(j) < 100
                count = count + 1;
                file_time=[rep,'/run_', num2str(run(j)), '/for_analysis/timevar'];
                full_file=importdata(file_time);
                timevar=full_file.data;
                t  = timevar(:,1);    
                em = timevar(:,4);    % dim + em
                ev = timevar(:,5);    % dim + ev
                
                % dim - parameter F, related to Re
                if ii==1
                    F = sqrt(3.0) * nuinv(j);
                    A = nuinv(j) ./ F;
                    C = sqrt(3.0) * nuinv(j) ./ F;
                else
                    F = sqrt( 1.0 + 1.25/nuinv(j)^2 ) * nuinv(j)^2;
                    A = 2.55 * nuinv(j) ./ F;
                    C = sqrt(0.125) * nuinv(j) ./ F;
                end
                % dim - d forcing lengthscale 
                % F,A,C include k = 2*pi/d=1
                % however, all the other quantities include d, not k
                d = 2 * pi;
                em = em.* nuinv(j)^2 / F * d^2;       % dim - em   
                ev = ev.* nuinv(j)^2 / F^2 * d^2;     % dim - ev  
                g = 1./sqrt(ev);                % dim - g = 1/sqrt(ev)
                k2 = em./ev;                  % dim - k^2 = em/ev
                Pm = etainv(j)/nuinv(j);   

                % for the approximation, create an array for the error bars
                lower = (g.^2 - C * g) / (F*A) -1.0;
                upper = (Pm .* (g - 1.0));
limx = [0 100];
                hFig = figure(count);
                set(hFig, 'Position', [100, 50, 600, 500]); 
                set(gca, 'FontSize', 14)
                    subplot(2,1,1)
                        p1= semilogy(t,em, ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','E_m'); 
                        title(['P_m=', num2str(Pm),'in [', num2str(box(jj)), ',2,1]']);
                        %xlim(limx)
                        %ylim([min(em(end),ev(end))*0.8 max(em(end),ev(end))*1.2])
                        hold on;
                        p2= semilogy(t,ev, ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',1.5,...
                            'DisplayName','E_k');
                        legend([p1 p2])%,'Location','west', 'Orientation','Horizontal')
                        
                    subplot(2,1,2)
                        p3 = plot(t,real(k2), ...
                            'LineStyle', '-',...
                            'color',color{1},...
                            'LineWidth',1.5);
                        xlabel('Time')
                        ylabel('E_m/E_k')
                        hold on;
                        p3 = plot(t,imag(k2), ...
                            'LineStyle', '--',...
                            'color',color{1},...
                            'LineWidth',1.5);
                        p3 = plot(t,real(lower), ...
                            'LineStyle', '-',...
                            'color',color{2},...
                            'LineWidth',1.5);
                        p3 = plot(t,imag(lower), ...
                            'LineStyle', '--',...
                            'color',color{2},...
                            'LineWidth',1.5);
                        p3 = plot(t,real(upper), ...
                            'LineStyle', '-',...
                            'color',color{3},...
                            'LineWidth',1.5);
                        p3 = plot(t,imag(upper), ...
                            'LineStyle', '--',...
                            'color',color{3},...
                            'LineWidth',1.5);
                        %xlim(limx)
            end
        end
    end
end