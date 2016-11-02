% First used for the mean-field approximation of u_III,
% which has alpha a(x) and hence, we set B=B(x,t).
% The resulting eigen matrix gives:
% eigenvalue: growth rate 
% eigenfunction: iBy + Bz = p = sum_n pn exp(inmx), 
% where m=1/Lx and a(x)=a0 cos^2(mx)
% see log_book_2



close all;
clear all;
% Provide values for constants involved:
% velocity profile: 
zi = sqrt(-1);
D_const=1;
ku=1;  
%arr_uii=[128 80 32 20 16 10 8];  
%rmc_uii=[1.530000e-01 2.350000e-01 0.3247 0.39 0.4647 0.5317 0.6225];
arr_uii=[8 10 16 20 32]; 
rmc_uii=[0.627 0.54 0.4647 0.35 0.27];

rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii/results/box_';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

files2=[{[rep1,'8/result_summary']}... 
        {[rep1,'10/result_summary']}...
        {[rep1,'16/result_summary']}...
        {[rep1,'20/result_summary']}...
        {[rep1,'32/result_summary']}];
    
legendInfo2 = [     {'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'} ];      
% matrix dimensions @n = [-@n_max, @n_max]
% @n is either positive or negative and so is n_max 
    % For boxes 32,20,16, rates dont change between 170~324
    % For boxes 32,20,8,2: 7-705 same 

nn=[11 16];
% box sizes to consider
ar = [8 10 16 20 32];
for ind=1:5
    color{ind}=[0,0,0];
end
% color{1} = [0.83,0,0.17]; 
% color{2} = [1,0.67,0.33]; 
% color{3} = [0.33,0.33,1.0];   
% color{4} = [0,0.5,0]; 
% color{5} = [0.5,0,0.5]; 
for co=1:2 
 
    n_max=nn(co);
    for box = 1:size(ar,2) 
        % color{box}=rand(1,3); 
        % wavenumber corresponding to the longest wavelength on the box
        m =1./ar(box);

        count=1; % data output count for each rm
        % magnetic diffusivity & magnetic Reynolds number
        for rm=0.001:0.001:1.0
        %for rm=10:1.0:500  
            eta = (1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))/rm;
            %eta = (D_const*sqrt(m^2/3.+0.5))/rm;

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
            [V,DD] = eig(D,'vector');
            [gr,Ig] = max(DD);
            [kk,Ik] = max(V(:,Ig));
            grr(count) = gr;
            rmm(count) = rm;
            kkk(count) = abs(n(Ik));
            count = count + 1;
        end
   
       
      
    figure(1)
        if (mod(co,2)~=0)
            subplot(2,size(ar,2) ,box)
                plot(rmm(grr>0),1./(kkk(grr>0)*m),'o',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.5);
                hold on;
                
                title(legendInfo2{box},'fontsize',16, 'Interpreter', 'latex');
                xlabel('R_m');
                if (box==1)
                    ylabel('L_{mean-field}/2\pi [L]');
                end
                 set(gca, 'FontSize', 14)     
           
            subplot(2,size(ar,2),box+size(ar,2))
                plot(rmm(grr>0),grr(grr>0),'o',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.5);
                hold on;
                  xlabel('R_m');
                  if (box==1)
                    ylabel('\sigma [T^{-1}_{turnover}]');
                end
                set(gca, 'FontSize', 14)     
                
            end

        if (mod(co,2)==0)
            subplot(2,size(ar,2) ,box)
                plot(rmm(grr>0),1./(kkk(grr>0)*m),'.',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.0);
                hold on;
                %title(legendInfo{box});
                xlabel('R_m');
                if (box==1)
                    ylabel('L_{mean field}/2\pi [L]');
                end
                 set(gca, 'FontSize', 14)
               
            subplot(2,size(ar,2),box+size(ar,2))
                plot(rmm(grr>0),grr(grr>0),'.',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.0);
                hold on;  
                xlabel('R_m');
                if (box==1)
                    ylabel('\sigma [T^{-1}_{turnover}]');
                end
                set(gca, 'FontSize', 14)     
            end
        end
    %legend(legendInfo,'Orientation','horizontal','Location','northoutside','fontsize',16, 'Interpreter', 'latex', 'Box', 'off');
    for j=1:size(files2,2)
        full_file=importdata(files2{1,j});
        timevar=full_file.data;
        nvar=size(timevar,2);
        nentry=n;
        n=size(timevar,1);
        % Create a table with all variables
        %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
        tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,2);
        Lx = tblB{:,1};
        etainv = tblB{:,2};
        Gr = tblB{:,3};
        lBx= tblB{:,4};
        lBy= tblB{:,5};
        lBz= tblB{:,6};
        
        m = 1./Lx(1); 
        Rm = (1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))*etainv;
        
       figure(1)
        subplot(2,size(files2,2) ,j)
            plot(Rm(Gr > -0.01),max(lBy(Gr > -0.01),lBz(Gr > -0.01))./(2*pi),'+',...
                    'LineStyle', 'none',...
                    'color','r',...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', 'r' , ...
                    'MarkerFaceColor','r', ...
                    'MarkerSize',5.0);
            hold on;
            plot([rmc_uii(j)*(1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2)) rmc_uii(j)*(1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))],[-0.01 max(ar)],'k:');
           xlim([0 rmm(end)]);
           ylim([-0.01 max(ar)]);
            
        
       subplot(2,size(files2,2) ,j+size(files2,2))
            plot(Rm(Gr > -0.01),Gr(Gr > -0.01),'+',...
                    'LineStyle', 'none',...
                    'color','r',...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor', 'r' , ...
                    'MarkerFaceColor', 'r', ...
                    'MarkerSize',5.0);
            hold on;
            plot([rmc_uii(j)*(1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2)) rmc_uii(j)*(1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2))],[-0.01 max(grr)],'k:');
            xlim([0 rmm(end)]);
            ylim([-0.01 max(grr)]);
    end
end