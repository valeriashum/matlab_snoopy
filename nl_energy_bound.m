close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep1=[{'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/box_', '/store/ASTRO/vs391/nonlinear_dynamo/u_iii/box_'}];
rep2=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc'     , '/store/ASTRO/vs391/kinematic_dynamo/u_iii'}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box=[10 16 32];
color{1}=[0.83 0.17 0.17];
color{2}=[0 0 0];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(rep1,2)
    clear tblA tblB run cas 
    for jj=1:size(box,2)
        clear em ev yneg ypos upper lower Pm k2 Ek nuinv etainv  
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
        
         rows = 4;
            if mod(size(nuinv,1),rows)==0
                cols = fix(size(nuinv,1)/rows);
            else 
                cols = fix(size(nuinv,1)/rows) + 1;
            end
        count = 0;
        for j=1:size(nuinv,1)
            clear em ev
            if cas(j) ~= 0 && run(j) < 100
                count = count + 1;
                file_time=[rep,'/run_', num2str(run(j)), '/for_analysis/timevar'];
                full_file=importdata(file_time);
                timevar=full_file.data;
                em = timevar(end,4);    % dim + em
                ev = timevar(end,5);    % dim + ev
                
                % dim - parameter F, related to Re
                if ii==1
                    F = sqrt(3) * nuinv(j);
                    A = nuinv(j) ./ F;
                    C = sqrt(3) * nuinv(j) ./ F;
                else
                    F = sqrt( 1.0 + 1.25/nuinv(j)^2 ) * nuinv(j)^2;
                    A = 2.55 * nuinv(j) ./ F;
                    C = sqrt(0.125) * nuinv(j) ./ F;
                end
                % dim - d forcing lengthscale 
                % F,A,C include k = 2*pi/d=1
                % however, all the other quantities include d, not k
                d = 1;%2 * pi;
                em = em * nuinv(j)^2 / F * d^2;       % dim - em   
                ev = ev * nuinv(j)^2 / F^2 * d^2;     % dim - ev  
                g = 1./sqrt(ev);                % dim - g = 1/sqrt(ev)
                k2(j) = real(em/ev);                  % dim - k^2 = em/ev
                Pm(j) = etainv(j)/nuinv(j);   
                Ek(j) = ev;
                % for the approximation, create an array for the error bars
                lower(j) = (g^2 - C * g) / (F*A) -1;
                yneg(j) = k2(j) - lower(j);
                upper(j) = Pm(j) * (g - 1);
                ypos(j) = k2(j) + upper(j);
                if  yneg(j) < 0
                    yneg(j) = NaN;
                end
                
            end
        end
        if count > 0
        hFig = figure(1);
        set(hFig, 'Position', [100, 50, 600, 500]); 
        set(gca, 'FontSize', 14)
            subplot(2,size(box,2), (ii-1)*3 + jj)
                color_er = [160 160 200]/255;
%                 p1 = errorbar(Pm,k2,yneg,ypos,'o', ...
%                     'LineStyle', 'none',...
%                     'color',color_er,...
%                     'LineWidth',1.5,...
%                     'MarkerEdgeColor', color_er, ...
%                     'MarkerFaceColor', color_er, ...
%                     'MarkerSize',2.5);
                %set(gca,'xscale','log','yscale','log')
                
                
                p2 =semilogx(Pm,k2,'o', ...
                    'LineStyle', 'none',...
                    'color',color{ii},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color{ii}, ...
                    'MarkerFaceColor', color{ii}, ...
                    'MarkerSize',4.5);
                hold on;
                if ii==1
                    title(['Box=[', num2str(box(jj)), ',2,1]'], 'Interpreter', 'latex');
                else
                    xlabel('$P_m$','Interpreter', 'latex');
                end
                if jj==1
                    ylabel('$E_m/E_k$','Interpreter', 'latex');
                end
                xlim([0 inf])
                ylim([min(min(lower, min(upper, k2))) max(max(lower,max( upper, k2)))]) 
                 p3 =plot(Pm,lower,'v', ...
                    'LineStyle', 'none',...
                    'color',color_er,...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color_er, ...
                    'MarkerFaceColor', color_er, ...
                    'MarkerSize',4.5);
                 p4 =plot(Pm,upper,'^', ...
                    'LineStyle', 'none',...
                    'color',color{ii},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color_er, ...
                    'MarkerFaceColor', color_er, ...
                    'MarkerSize',4.5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig = figure(2);
        set(hFig, 'Position', [100, 50, 600, 500]); 
        set(gca, 'FontSize', 14)
            subplot(2,size(box,2), (ii-1)*3 + jj)
                color_er = [160 160 200]/255;
%                 p1 = errorbar(Pm,k2,yneg,ypos,'o', ...
%                     'LineStyle', 'none',...
%                     'color',color_er,...
%                     'LineWidth',1.5,...
%                     'MarkerEdgeColor', color_er, ...
%                     'MarkerFaceColor', color_er, ...
%                     'MarkerSize',2.5);
                %set(gca,'xscale','log','yscale','log')
                
                
                p2 =semilogx(Ek,k2,'o', ...
                    'LineStyle', 'none',...
                    'color',color{ii},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color{ii}, ...
                    'MarkerFaceColor', color{ii}, ...
                    'MarkerSize',4.5);
                hold on;
                if ii==1
                    title(['Box=[', num2str(box(jj)), ',2,1]'], 'Interpreter', 'latex');
                else
                    xlabel('$E_k$','Interpreter', 'latex');
                end
                if jj==1
                    ylabel('$E_m/E_k$','Interpreter', 'latex');
                end
                xlim([0 inf])
                ylim([min(min(lower, min(upper, k2))) max(max(lower,max( upper, k2)))]) 
                 p3 =plot(Ek,lower,'v', ...
                    'LineStyle', 'none',...
                    'color',color_er,...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color_er, ...
                    'MarkerFaceColor', color_er, ...
                    'MarkerSize',4.5);
                 p4 =plot(Ek,upper,'^', ...
                    'LineStyle', 'none',...
                    'color',color{ii},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', color_er, ...
                    'MarkerFaceColor', color_er, ...
                    'MarkerSize',4.5);
        end
    end
end
 set(gca, 'FontSize', 14)