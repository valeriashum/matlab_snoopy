% KINEMATIC REGIME
% Plot the growth rates vs Rm

clear all;
close all;
rep1=[{'/store/ASTRO/vs391/kinematic_dynamo/u_abc'     , '/store/ASTRO/vs391/kinematic_dynamo/u_iii'}];
box=[8 10 16 20 32];


Rm_tr_abc = [12.949132605771180,13.290036428387761,5.715832712497174,4.209549700787671,4.075457990762875];
Rm_tr_iii = [2.400446928338508,2.721222454636513,2.478342661461882,2.196408715239213,2.061712447831615];
for ii=1:size(rep1,2)
    for jj=1:size(box,2)
        color{1}=[0.83 0.17 0.17];
        color{2}=[0 0 0];
        file_time=[rep1{ii},['/box_', num2str(box(jj)), '_modes']];
        full_file=importdata(file_time);
        timevar=full_file;
        
        % Create a table with all variables
        %tblA = table(Lx,Rm,Gr,kBx,kBy,kBz);
        tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,2);
        Lx = tblB{:,1};
        Rm = tblB{:,2};
        Gr = tblB{:,3};
        kx= tblB{:,4};
        ky= tblB{:,5};
        kz= tblB{:,6};        
    
        if ii == 2 
            m = 1.0/box(jj);
            Rm = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))) * Rm;
        end   
        
        % get rid of negative growth rates
        Lx = Lx(Gr > 0);
        Rm = Rm(Gr > 0);
        Gr = Gr(Gr > 0);
        kx = kx(Gr > 0);
        ky = ky(Gr > 0);
        kz = kz(Gr > 0);
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        hFig = figure(1);
        set(hFig, 'Position', [100, 50, 1100, 600]); 
        set(gca, 'FontSize', 14)
        subplot(2,size(box,2),(ii-1)*size(box,2) + jj)
            set(gca, 'FontSize', 16)
            loglog(Rm(Gr>0), Gr(Gr>0),'o',...
                        'LineStyle', 'none',...
                        'color',color{ii},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{ii}, ...
                        'MarkerFaceColor', color{ii}, ...
                        'MarkerSize',4.5); 
            xlim([0.1 50])  
            ylim([1e-5 1])
            hold on;
           if ii==1
               title(['Box=[', num2str(box(jj)),',2,1]'],'fontsize',16, 'Interpreter', 'latex')
           end
           if ii==2 
                xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
           end
           if jj==1
                ylabel('$\sigma [T^{-1}]$', 'Interpreter', 'latex')
           end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            if ii == 1
                Rm_tr = Rm_tr_abc;
            else
                Rm_tr = Rm_tr_iii;
            end
            color{2}=[0.83 0.17 0.17];
            color{1}=[0 0 0];
            %%%%%%%%%%%%%%%%
            p1 = polyfit(log(Rm(Rm < Rm_tr(jj))),log(Gr(Rm < Rm_tr(jj))),1); 
            m1= p1(1);
            b1 = exp(p1(2));
            loglog(Rm(Rm < Rm_tr(jj)), b1*Rm(Rm < Rm_tr(jj)).^m1, ...
                    'LineStyle', ':',...
                    'color',color{ii},...
                    'LineWidth',1.5);  
            str = ['$ \sim R_m^{ ',num2str(m1,'%1.1f'), '}$'];
            text(0.5,0.001,str, 'Interpreter', 'latex','fontsize',14)    
            %%%%%%%%%%%%%%%%
            p1 = polyfit(log(Rm(Rm > Rm_tr(jj))),log(Gr(Rm > Rm_tr(jj))),1); 
            m1= p1(1);
            b1 = exp(p1(2));
            loglog(Rm(Rm > Rm_tr(jj)), b1*Rm(Rm > Rm_tr(jj)).^m1, ...
                    'LineStyle', ':',...
                    'color',color{ii},...
                    'LineWidth',1.5); 
            str = ['$ \sim R_m^{ ',num2str(m1,'%1.1f'), '}$'];
            text(1.8,0.05,str, 'Interpreter', 'latex','fontsize',14)      
            %%%%%%%%%%%%%%%%
            
    end
    set(gca, 'FontSize', 14)
end