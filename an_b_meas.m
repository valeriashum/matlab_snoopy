close all;
clear all;

rep={'/store/ASTRO/vs391/kinematic_dynamo/u_iii','/store/ASTRO/vs391/kinematic_dynamo/u_abc'};
box=[8 10 16 20 32];
color{1} = [0.83,0,0.17]; 
color{2} = [0,0,0];

for ii = 1:size(rep,2)
    clear Rm_tr Rm2
    for jj = 1:size(box,2)
        file_1=[rep{ii},'/b_meas_', num2str(box(jj)), '.dat'];
        timevar1=importdata(file_1);
        tblA = table(timevar1(:,1),timevar1(:,2),timevar1(:,3));
        tblB = sortrows(tblA,1); 
        if ii == 1
            Rm = tblB{:,1};
        else
            Lx = box(jj);
            Rm = tblB{:,1}.*(1.5+1./Lx.^2)/sqrt(1.5*(3+2*1./Lx.^4+9*1./Lx.^2));
        end
        LB = tblB{:,2};
        Ratio = tblB{:,3};
 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig1 = figure(1);
        set(hFig1, 'Position', [100, 50, 1200, 345]); 
        set(gca, 'FontSize', 14)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(1,5,jj)
            loglog(Rm,LB,'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5);
            hold on;
            ylim([1 35*6.28])
            title(['Box=[', num2str(box(jj)), ',2,1]'], 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            if (jj==1)
                ylabel('Lengthscale [L]','fontsize',16, 'Interpreter', 'latex');
            end 
%         p = polyfit(log(Rm),log(LB),1); 
%         m1= p(1);
%         b1 = exp(p(2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig2 = figure(2);
        set(hFig2, 'Position', [100, 50, 1100, 900]); 
        set(gca, 'FontSize', 14)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,5,jj + 5*(ii-1))
            loglog(Rm,Ratio,'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5);     
            hold on;
            if ii == 1
            title(['Box=[', num2str(box(jj)), ',2,1]'], 'Interpreter', 'latex');
            end
            if (jj==1)
                ylabel('$ E_m^{large}/E_m^{small}$','fontsize',16, 'Interpreter', 'latex');
            end 
            plot([0.1 100],[1 1],'k:') 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        % Large Fit
                                                        Rm_l = Rm(Ratio > 1);
                                                        Ratio_l = Ratio(Ratio > 1);
                                                        p = polyfit(log(Rm_l),log(Ratio_l),1); 
                                                        m1= p(1);
                                                        b1 = exp(p(2));
%                                                                 loglog(Rm_l, b1*Rm_l.^m1, ...
%                                                                         'LineStyle', ':',...
%                                                                         'color',color{2},...
%                                                                         'LineWidth',1.5);
%                                                                 xlim([[0.1 100]])
%                                                                 str = ['$ \sim R_m^{ ',num2str(m1,'%1.1f'), '}$'];
%                                                                 text(1,100,str, 'Interpreter', 'latex','fontsize',14)
                                                        % Small Fit
                                                        Rm_s = Rm(Ratio < 1);
                                                        Ratio_s = Ratio(Ratio < 1);
                                                        p = polyfit(log(Rm_s),log(Ratio_s),1); 
                                                        m2 = p(1);
                                                        b2 = exp(p(2));
%                                                                 loglog(Rm_s, b2*Rm_s.^m2, ...
%                                                                         'LineStyle', ':',...
%                                                                         'color',color{2},...
%                                                                         'LineWidth',1.5);
%                                                                 xlim([[0.1 100]])
%                                                                 ylim([[0.1 10000]])
%                                                                 str = ['$ \sim R_m^{ ',num2str(m2,'%1.1f'), '}$'];
%                                                                 text(0.4,0.5,str, 'Interpreter', 'latex','fontsize',14)
                                                       % calculate Rm interesections
                                                        Rm1=(1./b1)^(1./m1);
                                                        Rm2=(1./b2)^(1./m2);
    
                                                        % Save to an array: 
                                                        Rm_tr(jj) = Rm2;
                                                        Rm_tr_abc(jj) = Rm_tr(jj);
                                                        %Rms=(b2/b1)^(1/(m1-m2));
                                                        e = std([Rm1 Rm2]);
                                                        %fprintf('For uABC in a box=%d,Rm_tr=%f\n',box(jj),Rm2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hFig3 = figure(3);
    set(hFig3, 'Position', [100, 50, 1100, 900]); 
    set(gca, 'FontSize', 14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(box(jj),Rm2,'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5);
        errorbar(box(jj),Rm2, e,'color',color{ii} )      
        xlabel('Aspect Ratio','fontsize',16, 'Interpreter', 'latex');
        if (jj==1)
            ylabel('$ Rm(E_m^l=E_m^s)$','fontsize',16, 'Interpreter', 'latex');
        end 
        hold on; 
        ratio1_abc(jj) = m1;
        ratio2_abc(jj) = m2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hFig4 = figure(4);
        set(hFig4, 'Position', [100, 50, 1200, 750]); 
        set(gca, 'FontSize', 14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,5,jj+(ii-1)*5)
            loglog(Rm,LB,'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5);
            hold on;
            if (ii==2)
            xlabel('$R_m$', 'Interpreter', 'latex');
            end
            if (jj==1)
                ylabel('Lengthscale [L]', 'Interpreter', 'latex');
            end 
            if ii == 1
                title(['Box=[', num2str(box(jj)), ',2,1]'], 'Interpreter', 'latex');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                                                            p1 = polyfit(log(Rm(Rm < Rm_tr(jj))),log(LB(Rm < Rm_tr(jj))),1); 
                                                            m1= p1(1);
                                                            b1 = exp(p1(2));
                                                                    loglog(Rm(Rm < Rm_tr(jj)), (b1+0.5)*Rm(Rm < Rm_tr(jj)).^(-1.0), ...
                                                                            'LineStyle', ':',...
                                                                            'color','k',...
                                                                            'LineWidth',1.5);
%                                                                     xlim([[0.1 50]])
%                                                                     str = ['$ \sim R_m^{ ',num2str(m1,'%1.1f'), '}$'];
%                                                                     text(1,20,str, 'Interpreter', 'latex','fontsize',14)    
                                                            p2 = polyfit(log(Rm(Rm > Rm_tr(jj))),log(LB(Rm > Rm_tr(jj))),1); 
                                                            m2= p2(1);
                                                            b2 = exp(p2(2));
                                                                    loglog(Rm(Rm > Rm_tr(jj)), (b2+0.5)*Rm(Rm > Rm_tr(jj)).^(-0.5), ...
                                                                            'LineStyle', ':',...
                                                                            'color','k',...
                                                                            'LineWidth',1.5);
%                                                                     xlim([[0.1 50]])
%                                                                     str = ['$ \sim R_m^{ ',num2str(m2,'%1.1f'), '}$'];
%                                                                     text(1,1.5,str, 'Interpreter', 'latex','fontsize',14)    
%                                                             m1_ar_iii(jj) = m1;
%                                                             m2_ar_iii(jj) = m2; 
                                                             set(gca, 'FontSize', 14)   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_2=[rep{ii},'/results/box_', num2str(box(jj)), '/result_summary'];       
    file2=importdata(file_2);
    timevar2=file2.data;
    tblA2 = table(timevar2(:,1),timevar2(:,2),timevar2(:,3),timevar2(:,4),timevar2(:,5),timevar2(:,6));
    % Sort the rows of the table based on Rm
    tblB2 = sortrows(tblA2,2);
    Lx = tblB2{:,1};
    if ii == 1
        Rm = tblB2{:,2};
    else
        Lx = box(jj);
        Rm = tblB2{:,2}.*(1.5+1./Lx.^2)/sqrt(1.5*(3+2*1./Lx.^4+9*1./Lx.^2));
    end
    Gr = tblB2{:,3};
    lBx= tblB2{:,4};
    lBy= tblB2{:,5};
    lBz= tblB2{:,6};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,5,jj)
        loglog(Rm,max(lBy,lBz),'^',...
                'LineStyle', 'none',...
                'color',color{ii},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{ii}, ...
                'MarkerFaceColor', color{ii}, ...
                'MarkerSize',4.5);
        hold on;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end            
end