
rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

files1=[%{[rep1,'/box_8_modes']}... 
        {[rep1,'/box_10_modes']}...
        {[rep1,'/box_16_modes']}...  
        {[rep1,'/box_20_modes']}... 
        {[rep1,'/box_32_modes']}];
        
files2=[%{[rep2,'/box_8_modes']}... 
        {[rep2,'/box_10_modes']}...
        {[rep2,'/box_16_modes']}...  
        {[rep2,'/box_20_modes']}... 
        {[rep2,'/box_32_modes']}];
    
legendInfo1 = [     %{'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'}];
D_const=1;
ku=1;
% color{1} = [0.83,0,0.17]; color{2} = [1,0.67,0.33]; color{3} = [0.33,0.33,1.0];   
% color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5]; 
for ind=1:5
    color{ind}=[0.83,0,0.17];
end

%%%%%%%%%
%%%%%%%%%
% U_III %
%%%%%%%%%
%%%%%%%%%

 for j=1:size(files1,2)
    
    full_file=importdata(files1{1,j});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1)
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,2);
    Lx = tblB{:,1};
    Rm = tblB{:,2}.*(1.5+1./Lx(1).^2)/sqrt(1.5*(1+2*1./Lx(1).^4+9*1./Lx(1).^2));
    Gr = tblB{:,3};
    kx= tblB{:,4};
    ky= tblB{:,5};
    kz= tblB{:,6};
    
    % save critical Rm into small scales
    
    for ind=1:ndat
        if (ky(ind) == 0)
            if (kx(ind) <= 0.25)
                rm_0p2(j) = Rm(ind);
                box_0p2(j) = Lx(1);
            end
            if (kx(ind) <= 0.1)
                rm_0p1(j) = Rm(ind);
                box_0p1(j) = Lx(1);
            end
            if (kx(ind) <= 0.4)
                rm_0p4(j) = Rm(ind);
                box_0p4(j) = Lx(1);
            end
        end
         % dynamo onset
         Rm_pos=Rm(Gr>0);
         Gr_pos=Gr(Gr>0);
         rmc(j)= Rm_pos(1);
         box_rmc(j) = Lx(1);
    end
       
    figure(3)
    subplot(1,size(files1,2),j)
    set(gca, 'FontSize', 18)    
            plot(Rm(Gr>0),1./kx(Gr>0),'o',...
                'LineStyle', 'none',...
                'color','k',...
                'LineWidth',1.5,...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'k', ...
                'MarkerSize',6.5);
            hold on;
            title(legendInfo1{j},'fontsize',16, 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            if (j==1)
                ylabel('$1/k_x [L^{-1}]$','fontsize',16, 'Interpreter', 'latex');
            end
            xlim([0 40]);
            ylim([0 32]);  
    
    figure(1)
        subplot(2,size(files1,2),j)
            plot(Rm(ky==0 & Gr>0),1./kx(ky==0 & Gr>0),'o',...
                'LineStyle', 'none',...
                'color',color{j},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',6.5);
            hold on;
            title(legendInfo1{j},'fontsize',16, 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            if (j==1)
                ylabel('$1/k_x [L^{-1}]$','fontsize',16, 'Interpreter', 'latex');
            end
            xlim([0 2.5]);
            ylim([0 32]);  
        
            
 end  
 colory{1}=[0.45 0.45 0.45];
 colory{2}=[0.91 0.91 0.87];
 colory{3}=[0.61 0.61 0.68];
 
 set(gca, 'FontSize', 16)
 figure(4)
         %subplot(2,size(files1,2), size(files1,2)+1:size(files1,2)*2)
         subplot(1,2,2) 
         plot(box_0p4(rm_0p4>0), rm_0p4(rm_0p4>0),'o',...
             'LineStyle', ':',...
                'color',colory{1},...
                'LineWidth',1.5,...
                'MarkerEdgeColor',colory{1} , ...
                'MarkerFaceColor', colory{1}, ...
                'MarkerSize',2.5); 
            xlabel('$L_x/2\pi$','fontsize',16, 'Interpreter', 'latex'); 
            ylabel('$R_{m}^{trans}(k_x)$','fontsize',16, 'Interpreter', 'latex')
            title('$u_{III}$','fontsize',16, 'Interpreter', 'latex')
            xlim([10 32])
            ylim([0 2.5])
            hold on;
            area(box_0p4(rm_0p4>0), rm_0p4(rm_0p4>0), 'FaceColor' ,colory{1});
          set(gca, 'FontSize', 20)
            
            plot(box_0p2(rm_0p2>0), rm_0p2(rm_0p2>0),'o',...
                'LineStyle', ':',...
                'color',colory{2},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', colory{2}, ...
                'MarkerFaceColor', colory{2}, ...
                'MarkerSize',2.5);
            area(box_0p2(rm_0p2>0), rm_0p2(rm_0p2>0),'FaceColor', colory{2});
            set(gca, 'FontSize', 20)
            plot(box_0p1(rm_0p1>0), rm_0p1(rm_0p1>0),'o',...
                'LineStyle', ':',...
                'color',colory{3},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', colory{3}, ...
                'MarkerFaceColor', colory{3}, ...
                'MarkerSize',2.5);
            area(box_0p1(rm_0p1>0), rm_0p1(rm_0p1>0),'FaceColor' ,colory{3});
            area(box_rmc, rmc,'FaceColor', 'w');
             text(27, 0.55, '$k_x < 0.1$','color','k','fontsize',16, 'Interpreter', 'latex' )
             text(27, 1.62, '$k_x < 0.25$','color','k','fontsize',16 , 'Interpreter', 'latex')
             text(10, 1.62, '$k_x < 0.4$','color','k','fontsize',16, 'Interpreter', 'latex' )
             set(gca, 'FontSize', 20)
% approximation              
nn=[11 16];
% box sizes to consider
ar = [8 10 16 20 32];
for ind=1:5
    color{ind}=[0,0,0];
end

for co=1:2 
    n_max=nn(co);
    for box = 1:size(ar,2) 
        % color{box}=rand(1,3); 
        % wavenumber corresponding to the longest wavelength on the box
        m =1./ar(box);

        count=1; % data output count for each rm
        % magnetic diffusivity & magnetic Reynolds number
        for rm=0.001:0.001:3.0
        %for rm=10:1.0:500  
            eta = (1.5+m^2)/sqrt(1.5*(1+2*m^4+9*m^2))/rm;
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
            subplot(2,size(files1,2),box)
                plot(rmm(grr>0),1./(kkk(grr>0)*m),'o',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.5);
                hold on;    set(gca, 'FontSize', 20) 
        end

        if (mod(co,2)==0)
            subplot(2,size(files1,2),box)
                plot(rmm(grr>0),1./(kkk(grr>0)*m),'.',...
                        'LineStyle', 'none',...
                        'color',color{box},...
                        'LineWidth',1.5,...
                        'MarkerEdgeColor', color{box} , ...
                        'MarkerFaceColor', color{box}, ...
                        'MarkerSize',2.0);
                hold on;    set(gca, 'FontSize', 20)
            end
        end                  
end


%%%%%%%%%
%%%%%%%%%
% U_ABC %
%%%%%%%%%
%%%%%%%%%
for ind=1:5
    color{ind}=[0.83,0,0.17];
end
 for j=1:size(files2,2)
    
    full_file=importdata(files2{1,j});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,2);
    Lx = tblB{:,1};
    Rm = tblB{:,2};
    Gr = tblB{:,3};
    kx= tblB{:,4};
    ky= tblB{:,5};
    kz= tblB{:,6};
    
    % save critical Rm into small scales
    
    for ind=1:ndat
        if (ky(ind) == 0)
            if (kx(ind) <= 0.25)
                rm_0p2(j) = Rm(ind);
                box_0p2(j) = Lx(1);
            end
            if (kx(ind) <= 0.1)
                rm_0p1(j) = Rm(ind);
                box_0p1(j) = Lx(1);
            end
            if (kx(ind) <= 0.4)
                rm_0p4(j) = Rm(ind);
                box_0p4(j) = Lx(1);
            end
        end
         % dynamo onset
         Rm_pos=Rm(Gr>0);
         Gr_pos=Gr(Gr>0);
         rmc(j)= Rm_pos(1);
         box_rmc(j) = Lx(1);
    end
    figure(3)
        subplot(1,size(files2,2),j)
            plot(Rm(Gr>0),1./kx(Gr>0),'o',...
                'LineStyle', 'none',...
                'color','r',...
                'LineWidth',1.5,...
                'MarkerEdgeColor', 'r', ...
                'MarkerFaceColor', 'r', ...
                'MarkerSize',6.5);
            hold on;
            plot([0 40],[2 2],'k:')
            plot([0 40],[1./(0.5+1/Lx(1)) 1./(0.5+1/Lx(1))],'k:')
            plot([0 40],[1./(0.5-1/Lx(1)) 1./(0.5-1/Lx(1))],'k:')
    figure(2)
    set(gca, 'FontSize', 20)
        subplot(2,size(files2,2),j)
            plot(Rm(ky==0 & Gr>0),1./kx(ky==0 & Gr>0),'o',...
                'LineStyle', 'none',...
                'color',color{j},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', color{j}, ...
                'MarkerFaceColor', color{j}, ...
                'MarkerSize',6.5);
            hold on;
            title(legendInfo1{j},'fontsize',16, 'Interpreter', 'latex');
            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
            if (j==1)
                ylabel('$1/k_x [L^{-1}]$','fontsize',16, 'Interpreter', 'latex');
            end
            xlim([0 2.5]);
            ylim([0 32]);  
        
            
 end  
 colory{1}=[0.45 0.45 0.45];
 colory{2}=[0.91 0.91 0.87];
 colory{3}=[0.61 0.61 0.68];
 figure(4)
 set(gca, 'FontSize', 16)
         %subplot(2,size(files2,2), size(files2,2)+1:size(files2,2)*2)
         subplot(1,2,1) 
         plot(box_0p4(rm_0p4>0), rm_0p4(rm_0p4>0),'o',...
             'LineStyle', ':',...
                'color',colory{1},...
                'LineWidth',1.5,...
                'MarkerEdgeColor',colory{1} , ...
                'MarkerFaceColor', colory{1}, ...
                'MarkerSize',2.5); 
            xlabel('$L_x/2\pi$','fontsize',16, 'Interpreter', 'latex'); 
            ylabel('$R_{m}^{trans}(k_x)$','fontsize',16, 'Interpreter', 'latex')
            title('$u_{abc}$','fontsize',16, 'Interpreter', 'latex')
            xlim([10 32])
            ylim([0 2.5])
            hold on;
            area(box_0p4(rm_0p4>0), rm_0p4(rm_0p4>0), 'FaceColor' ,colory{1});
          set(gca, 'FontSize', 20)
            
            plot(box_0p2(rm_0p2>0), rm_0p2(rm_0p2>0),'o',...
                'LineStyle', ':',...
                'color',colory{2},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', colory{2}, ...
                'MarkerFaceColor', colory{2}, ...
                'MarkerSize',2.5);
            area(box_0p2(rm_0p2>0), rm_0p2(rm_0p2>0),'FaceColor', colory{2});
            set(gca, 'FontSize', 20)
           
            plot(box_0p1(rm_0p1>0), rm_0p1(rm_0p1>0),'o',...
                'LineStyle', ':',...
                'color',colory{3},...
                'LineWidth',1.5,...
                'MarkerEdgeColor', colory{3}, ...
                'MarkerFaceColor', colory{3}, ...
                'MarkerSize',2.5);
            area(box_0p1(rm_0p1>0), rm_0p1(rm_0p1>0),'FaceColor' ,colory{3});
            area(box_rmc, rmc,'FaceColor', 'w');
             text(27, 0.55, '$k_x < 0.1$','color','k','fontsize',16, 'Interpreter', 'latex' )
             text(27, 1.62, '$k_x < 0.25$','color','k','fontsize',16 , 'Interpreter', 'latex')
             text(10, 1.62, '$k_x < 0.4$','color','k','fontsize',16, 'Interpreter', 'latex' )
             set(gca, 'FontSize', 20)


