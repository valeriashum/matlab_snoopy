cd % First used for the mean-field approximation of u_III,
% which has alpha a(x) and hence, we set B=B(x,t).
% The resulting eigen matrix gives:
% eigenvalue: growth rate 
% eigenfunction: iBy + Bz = p = sum_n pn exp(inmx), 
% where m=1/Lx and a(x)=a0 cos^2(mx)
% see log_book_2
clear all;
%close all;
rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';

files1=[{[rep1,'/box_32_modes']}... 
        {[rep1,'/box_20_modes']}...
        {[rep1,'/box_16_modes']}...  
        {[rep1,'/box_10_modes']}... 
        {[rep1,'/box_8_modes']}];
        



% Provide values for constants involved:
% velocity profile: 
D_const=1;
ku=1;  
% abc case in boxes of x:1:1 for smaller boxes
arr_111=[5 4 2 1];
rmc_111=[0.7805 0.9785 4.425 8.9];
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

nn=[111 210];
% box sizes to consider
ar = [32 20 16 10 8];   
%ar = [1000 100 10]
color{1} = [0.83,0,0.17]; color{2} = [1,0.67,0.33]; color{3} = [0.33,0.33,1.0];   
color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5]; color{6} = [0.5,0,0.5]; color{7} = [0.5,0,0.5]; color{8} = [0.5,0,0.5]; 

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
        for rm=0.0001:0.001:5.0
            eta = ((1.5+m^2)/sqrt(1.5*(3+2*m^4+9*m^2)))/rm;
           
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
            [V,DD] = eig(-zi*D,'vector');   %DD is a diagonal matrix with eigenvalues 
            V = real(V);                    %with a corresponding eigenvector V (column)
            [gr,Ig] = max(DD);
            [kk,Ik] = max(V(:,Ig));
            rgr_max(count) = gr;
            rmm(count) = rm;
            kkk(count) = abs(n(Ik));
            rm_entry(count)= rm;        
            
            count=count+1;
            
        end   % end rm loop 
        
        legendInfo{count_box} = ['box = [' num2str(ar(box)) ',2,1]'];
        fprintf('m=%d\n', ar(box));
        
        % Critical Rmc 
            [min_gr, I] = min(rgr_max(rgr_max>0));
            rm_entry_positive = rm_entry(rgr_max>0);
            rmc(box) = rm_entry_positive(I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        hFig = figure(1);
        set(hFig, 'Position', [100, 60, 800, 500]);
        set(gca, 'FontSize', 14)
        if (mod(n_max,2)~=0)   %odd or even modes
            %plot Rm vs max(rgr) and max(rgr)
            subplot(2,2,1)
            plot(rm_entry, rgr_max,'o',...
                        'LineStyle', 'none',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',3.0); 
            
            ylabel('\sigma_{im}'); 
            %ylabel('\sigma [L^2T^{-1}_{turnover}] '); 
            xlim([0 rm_entry(end)])
            %ylim([0 0.02]);
            hold on; 
            
            title({[' ODD modes [-',num2str(n_max),',',num2str(n_max),']']},'interpreter','Latex');
            subplot(2,2,3)
            plot(rm_entry,1./(kkk*m),'o',...
                        'LineStyle', 'none',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',3.0); 
            
                hold on;
            xlabel('R_m');
            ylabel('L');
        else

            subplot(2,2,2)
            set(gca, 'FontSize', 18)
                semilogy(rm_entry, rgr_max,'o',...
                        'LineStyle', '-',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',1.0);  
             set(gca, 'FontSize', 14)        
            %xlabel('R_m');
            %ylabel('\sigma_{im}'); 
            xlim([0 rm_entry(end)])
            %ylim([0 0.02]);
            hold on; 
            title({['EVEN modes [2,',num2str(n_max),']']},'interpreter','Latex');  
            legend(legendInfo,'fontsize',16,'Location','southeast', 'Interpreter', 'latex','Box','off')
            subplot(2,2,4)
                plot(rm_entry,1./(kkk*m),'o',...
                        'LineStyle', 'none',...
                        'color',color{count_box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{count_box}, ...
                        'MarkerFaceColor', color{count_box}, ...
                        'MarkerSize',3.0); 
            
                hold on;
            xlabel('R_m');
            ylabel('L');
        end
        count_box = count_box + 1;
    end % end m loop 
end

% figure(1)
%     subplot(1,2,1)
%         plot([0 1], [0 0], 'k:')
%     subplot(1,2,2)
%         plot([0 1], [0 0], 'k:') 
       
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
    Rm = tblB{:,2}*(1.5+1./Lx(1).^2)/sqrt(1.5*(3+2*1./Lx(1).^4+9*1./Lx(1).^2));
    Gr = tblB{:,3};
    kx= tblB{:,4};
    ky= tblB{:,5};
    kz= tblB{:,6};        
    
%     figure(1)
%         subplot(2,2,1)
%         set(gca, 'FontSize', 14)
%             plot(Rm(Gr > 0), Gr(Gr > 0),'o',...
%                         'LineStyle', 'none',...
%                         'color',color{j},...
%                         'LineWidth',1.5,...
%                         'MarkerEdgeColor',color{j} , ...
%                         'MarkerFaceColor', color{j}, ...
%                         'MarkerSize',4.0); 
%            hold on;
       
end
        figure(1)
        suptitle('Steady oscillatory solution for \alpha (x)= cos^2(mx) and B(x)')
        