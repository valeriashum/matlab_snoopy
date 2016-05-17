% First used for the mean-field approximation of u_III,
% which has alpha a(x) and hence, we set B=B(x,t).
% The resulting eigen matrix gives:
% eigenvalue: growth rate 
% eigenfunction: iBy + Bz = p = sum_n pn exp(inmx), 
% where m=1/Lx and a(x)=a0 cos^2(mx)
% see log_book_2
clear all;

rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

files1=[{[rep1,'/box_32_modes']}... 
        {[rep1,'/box_20_modes']}...
        {[rep1,'/box_16_modes']}...  
        {[rep1,'/box_10_modes']}... 
        {[rep1,'/box_8_modes']}];
        
files2=[{[rep2,'/box_32_modes']}... 
        {[rep2,'/box_20_modes']}...
        {[rep2,'/box_16_modes']}...  
        {[rep2,'/box_10_modes']}... 
        {[rep2,'/box_8_modes']}];



% Provide values for constants involved:
% velocity profile: 
D_const=1;
ku=1;  
%arr_uii=[128 80 32 20 16 10 8];  
%rmc_uii=[1.530000e-01 2.350000e-01 0.3247 0.39 0.4647 0.5317 0.6225];
arr_uii=[32 20 16 10 8]; 
rmc_uii=[0.27 0.35 0.4 0.54 0.627]; %REAL
%rmc_uii=[0.27 0.34 0.36 0.45 0.5];   % AIM
arr_211=[32 20 16 10 8];
rmc_211=[0.1753  0.2354   0.2695    0.3385   4.275000e-01];

% matrix dimensions @n = [-@n_max, @n_max]
% @n is either positive or negative and so is n_max 
    % For boxes 32,20,16, rates dont change between 170~324
    % For boxes 32,20,8,2: 7-705 same 

nn=[7 16];
% box sizes to consider
ar = [32 20 16 10 8];    
color{1} = [0.83,0,0.17]; color{2} = [1,0.67,0.33]; color{3} = [0.33,0.33,1.0];   
color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5]; color{6} = [0.5,0,0.5]; 

for co=1:2 
    n_max=nn(co);
    count_box=1;
    rmc = zeros(1,size(ar,2));
    for box = 1:size(ar,2) 
        %color{box} = rand(1,3);
        % wavenumber corresponding to the longest wavelength on the box
        m =1./ar(box);

        count=1; % data output count for each rm
        % magnetic diffusivity & magnetic Reynolds number
        for rm=0.0001:0.001:1.0
            eta = ((1.5+m^2)/sqrt(1.5*(1+2*m^4+9*m^2)))/rm;
           
            % @a0 = @D^2@ku/@eta
            a0=D_const^2*ku/eta;

            % Set imaginary i=zi
            zi=sqrt(-1);

            %Create diagonal matrix
            k_m=0;  % count for main diagonal 
            k_u=0;  % count for upper diagonal 
            k_l=0;  % count for lower diagonal


            % Set up the matrix 
            if (mod(n_max,2)~=0)
                % Preallocate arrays for diagonals 
                diag_main = zeros(1,n_max+1);
                diag_upper= zeros(1,n_max);
                diag_lower= zeros(1,n_max);

                for i=1:2:n_max*2+1
                    k_m=k_m+1;
                    n(k_m)=i-n_max-1;
                    %fprintf('Rm=%f, k_m=%d, n(k_m)=%f\n',rm, k_m,n(k_m))
                    diag_main(k_m)= a0*n(k_m)*m/2. - eta*n(k_m)^2*m^2;

                    if (i>1)
                        k_u=k_u+1;
                        diag_upper(k_u) = a0*(n(k_m)*m/4. - m/2.); 
                    end
                    if (i< n_max*2+1)
                        k_l=k_l+1;
                        diag_lower(k_l) = a0*(n(k_m)*m/4 + m/2.);
                    end

                    D_main = diag(diag_main); 
                    D_upper = diag(diag_upper, 1);
                    D_lower = diag(diag_lower, -1);
                    D=D_main +D_upper + D_lower;
                end

            else
                % Preallocate arrays for diagonals 
                diag_main = zeros(1,n_max/2);
                diag_upper= zeros(1,n_max/2-1);
                diag_lower= zeros(1,n_max/2-1);
                for i=2:2:n_max
                    k_m=k_m+1;
                    n(k_m)=i;
                    %fprintf('Rm=%f, k_m=%d, n(k_m)=%f\n',rm, k_m,n(k_m))
                    diag_main(k_m)=  a0*n(k_m)*m/2. - eta*n(k_m)^2*m^2;

                    if (i>2)
                        k_u=k_u+1;
                        diag_upper(k_u) = a0*(n(k_m)*m/4. - m/2.);
                    end
                    if (i< n_max)
                        k_l=k_l+1;
                        diag_lower(k_l) = a0*(n(k_m)*m/4 + m/2.);
                    end

                    D_main = diag(diag_main); 
                    D_upper = diag(diag_upper, 1);
                    D_lower = diag(diag_lower, -1);
                    D=D_main +D_upper + D_lower;  

                end
            end

            % Find the eigenvalues of D
            e = eig(D);


            % save data: 
            rm_entry(count)= rm;
            [rgr_max(count),I] = max(real(e));
            cgr_max(count) = imag(e(I));
            
            
            rm_abc(count)=1.0/eta;
            p_abc(count)=  m/rm_abc(count) -m^2/rm_abc(count)^2;
            
            
            count=count+1;
            
        end   % end rm loop 
        
        legendInfo{count_box} = ['box = [' num2str(ar(box)) ',2,1]'];
        fprintf('m=%d\n', ar(box));
        
        % Critical Rmc 
            [min_gr, I] = min(rgr_max(rgr_max>0));
            rm_entry_positive = rm_entry(rgr_max>0);
            rmc(box) = rm_entry_positive(I);
      figure(1)
      set(gca, 'FontSize', 18)
        if (mod(n_max,2)~=0)
            %plot Rm vs max(rgr) and max(rgr)
            subplot(2,2,1)
            set(gca, 'FontSize', 14)
            plot(rm_entry, rgr_max,'o',...
                        'LineStyle', 'none',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',1.0); 
            xlabel('R_m');
            ylabel('\sigma [L^2T^{-1}_{turnover}] '); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.02]);
            hold on; 
            legend(legendInfo,'Location','southeast','fontsize',16, 'Interpreter', 'latex','Box','off')
            title({['Asymptotic Approximation $\alpha(x)$ and $B(x)$ for ODD modes [-',num2str(n_max),',',num2str(n_max),']']},'interpreter','Latex', 'fontsize',18);
            
%             subplot(3,2,3)
%             plot(rm_entry, cgr_max, 'color',color{count_box},'LineWidth',1.5)
%             xlabel('1/\eta');
%             ylabel('Complex \sigma');
%             xlim([0 rm_entry(end)])
%             hold on;
%             count_box=count_box+1;
%             %legend(legendInfo,'Location','northwest')
            
            subplot(2,2,3:4)
            set(gca, 'FontSize', 18)
            ar_odd = ar(rmc>0);
            rmc_odd= rmc(rmc>0);
            p1=plot(ar(rmc>0), rmc(rmc>0),'o',...
                        'LineStyle', ':',...
                        'color','k',...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceColor', 'k', ...
                        'MarkerSize',4.5,...
                        'DisplayName','odd modes');
            ylabel('Critical R_m')
            xlabel('Aspect ratio L_x/L_z')
            hold on;
             set(gca, 'FontSize', 14)
        else

            subplot(2,2,2)
            set(gca, 'FontSize', 18)
            plot(rm_entry, rgr_max,'o',...
                        'LineStyle', '-',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',1.0);  
             set(gca, 'FontSize', 14)        
            xlabel('R_m');
            ylabel('\sigma [L^2T^{-1}_{turnover}]'); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.02]);
            hold on; 
            %legend(legendInfo,'Location','outsidesouth')
            title({['Asymptotic Approximation  $\alpha(x)$ and $B(x)$ for EVEN modes [2,',num2str(n_max),']']},'interpreter','Latex', 'fontsize',18);  

%             subplot(3,2,4)
%             plot(rm_entry, cgr_max, 'color',color{count_box},'LineWidth',1.5)
%             xlabel('1/\eta');
%             ylabel('Complex \sigma');
%             xlim([0 rm_entry(end)])
%             hold on;
             %count_box=count_box+1; 
%             %legend(legendInfo,'Location','northwest')
            
            subplot(2,2,3:4)
            set(gca, 'FontSize', 18)
            %area([0 max(ar)],[0.21 0.21], 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
            %text(1, 0.1, 'MEANFIELD', 'color', 'k')
            xlim([0 max(ar)])
            alpha(.02)
            hold on;
            ar_even = ar(rmc>0);
            rmc_even= rmc(rmc>0);
            p2=plot(ar(rmc>0), rmc(rmc>0),'o',...
                        'LineStyle', ':',...
                        'color','r',...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor','r', ...
                        'MarkerFaceColor', 'r', ...
                        'MarkerSize',4.5,...
                        'DisplayName','even modes');
            % From Snoopy, critical Rm :
            ar_snoopy= [10 8];
            rmc_snoopy=[1.1218 1.1458];
%             %p3=plot(ar_snoopy, rmc_snoopy,'o',...
%                         'LineStyle', ':',...
%                         'color','r',...
%                         'LineWidth',1.5,...
%                         'MarkerEdgeColor', 'r', ...
%                         'MarkerFaceColor', 'r', ...
%                         'MarkerSize',2.5,...
%                         'DisplayName','SNOOPY');
            legend([p1 p2])
            set(gca, 'FontSize', 14)
        end
        
        
         figure(4)
        set(gca, 'FontSize', 18)
        if (mod(n_max,2)~=0)
            %plot Rm vs max(rgr) and max(rgr)
            subplot(1,2,1)
            set(gca, 'FontSize', 18)
            plot(rm_entry, rgr_max,'o',...
                        'LineStyle', 'none',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',1.0); 
            xlabel('$R_m$', 'Interpreter', 'latex');
            ylabel('$\sigma$ [$L^2T^{-1}_{turnover}$] ', 'Interpreter', 'latex'); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.12]);
            hold on; 
            %legend(legendInfo,'Location','southeast','fontsize',16, 'Interpreter', 'latex','Box','off')
            title({['Asymptotic Approximation $\alpha(x)$ and $B(x)$ for ODD modes [-',num2str(n_max),',',num2str(n_max),']']},'interpreter','Latex', 'fontsize',18);
            
%              plot(rm_abc,p_abc, ...
%                 'LineStyle', ':',...
%                 'color',color{count_box}, ...
%                 'LineWidth',2.5);
            count_box=count_box+1;
        else

            subplot(1,2,2)
            set(gca, 'FontSize', 18)
            plot(rm_entry, rgr_max,'o',...
                        'LineStyle', '-',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',1.0);  
            xlabel('$R_m$', 'Interpreter', 'latex');
            ylabel('$\sigma$ $[L^2T^{-1}_{turnover}]$', 'Interpreter', 'latex'); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.12]);
            hold on; 
            legend(legendInfo,'Location','northeast','fontsize',20, 'Interpreter', 'latex','Box','off')
            %legend(legendInfo,'Location','outsidesouth')
            title({['Asymptotic Approximation  $\alpha(x)$ and $B(x)$ for EVEN modes [2,',num2str(n_max),']']},'interpreter','Latex', 'fontsize',18);  
            
            count_box=count_box+1; 
        end
    end % end m loop 
end

figure(2)
p1 = loglog(ar_odd, rmc_odd, ... 
    'LineStyle', '--',...
    'color','k',...
    'LineWidth',1.5,...
    'DisplayName','u_{iii}^{odd} approximation');

title('Critical $R_m$ in boxes of $[L_x,2,1]$','fontsize',18, 'Interpreter', 'latex');
set(gca, 'FontSize', 14)
xlabel('L_x [L]');
ylabel('R_{mc}');
hold on;
p2 = loglog(ar_even, rmc_even, ...
    'LineStyle', '-.',...
    'color','k',...
    'LineWidth',1.0,...
    'DisplayName','u_{iii}^{even} approximation');
p3 = loglog(ar, sqrt(1./ar), ...
    'LineStyle', '--',...
    'color','r',...
    'LineWidth',1.5,...
    'DisplayName','u_{abc} approximation');
 set(gca, 'FontSize', 14)
p4 = loglog(arr_uii, rmc_uii*(1.5+m^2)/sqrt(1.5*(1+2*m^4+9*m^2)), 'o',...
    'LineStyle', 'none',...
    'color','k',...
    'LineWidth',1.5,...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize',4.5,...
    'DisplayName','u_{iii}');
set(gca, 'FontSize', 14)
p5 = loglog(arr_211, rmc_211, 'o',...
    'LineStyle', 'none',...
    'color','r',...
    'LineWidth',1.5,...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize',4.5,...
    'DisplayName','u_{abc}');
set(gca, 'FontSize', 14)
legend([p1 p2 p3 p4 p5])
  
figure(1)
    subplot(2,2,1)
        plot([0 1], [0 0], 'k:')
    subplot(2,2,2)
        plot([0 1], [0 0], 'k:') 
  figure(4)
    subplot(1,2,1)
        plot([0 1], [0 0], 'k:')
    subplot(1,2,2)
        plot([0 1], [0 0], 'k:')          
 figure(3)
 plot(ar_even, rmc_odd.*sqrt(ar_even))
       
for j=1:size(files1,2)
    
    full_file=importdata(files1{1,j});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,2);
    Lx = tblB{:,1};
    Rm = tblB{:,2}*(1.5+1./Lx(1).^2)/sqrt(1.5*(1+2*1./Lx(1).^4+9*1./Lx(1).^2));
    Gr = tblB{:,3};
    kx= tblB{:,4};
    ky= tblB{:,5};
    kz= tblB{:,6};        
    
    figure(4)
    subplot(1,2,1)
            set(gca, 'FontSize', 18)
            plot(Rm, Gr,'o',...
                        'LineStyle', 'none',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor','k', ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.0); 
           hold on;
    subplot(1,2,2)
            set(gca, 'FontSize', 18)
            plot(Rm, Gr,'o',...
                        'LineStyle', 'none',...
                        'color',color{j},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceColor', color{j}, ...
                        'MarkerSize',4.0); 
           hold on;
end
        
        