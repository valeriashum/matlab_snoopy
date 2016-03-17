% First used for the mean-field approximation of u_III,
% which has alpha a(x) and hence, we set B=B(x,t).
% The resulting eigen matrix gives:
% eigenvalue: growth rate 
% eigenfunction: iBy + Bz = p = sum_n pn exp(inmx), 
% where m=1/Lx and a(x)=a0 cos^2(mx)
% see log_book_2




clear all;
% Provide values for constants involved:
% velocity profile: 
D_const=1;
ku=1;  

% matrix dimensions @n = [-@n_max, @n_max]
% @n is either positive or negative and so is n_max 
    % For boxes 32,20,16, rates dont change between 170~324
    % For boxes 32,20,8,2: 7-705 same 

nn=[7 16];
% box sizes to consider
ar = [128 80 32 20 16 10 8 5 4 2 1];    
for co1=1:size(ar,2)
    color{co1}=rand(1,3); 
end
for co=1:2 
    n_max=nn(co);
    count_box=1;
    rmc = zeros(1,size(ar,2));
    for box = 1:size(ar,2) 
        % wavenumber corresponding to the longest wavelength on the box
        m =1./ar(box);

        count=1; % data output count for each rm
        % magnetic diffusivity & magnetic Reynolds number
        for rm=0.0001:0.001:1.0
            eta = (D_const*sqrt(m^2+1.5))/rm;

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
            
            count=count+1;
        end   % end rm loop 
        
        legendInfo{count_box} = ['box = ' num2str(ar(box))];
        fprintf('m=%d\n', ar(box));
        
        % Critical Rmc 
            [min_gr, I] = min(abs(rgr_max));
            rmc(box) = rm_entry(I);
      
        if (mod(n_max,2)~=0)
            %plot Rm vs max(rgr) and max(rgr)
            subplot(3,2,1)
            plot(rm_entry, rgr_max, 'color',color{count_box},'LineWidth',1.5)
            xlabel('1/\eta');
            ylabel('Real \sigma'); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.02]);
            hold on; 
            legend(legendInfo,'Location','northeastoutside')
            title(['Asymptotic Approximation for a(x) and B(x) for ODD modes [-',num2str(n_max),',',num2str(n_max),']']);

            subplot(3,2,3)
            plot(rm_entry, cgr_max, 'color',color{count_box},'LineWidth',1.5)
            xlabel('1/\eta');
            ylabel('Complex \sigma');
            xlim([0 rm_entry(end)])
            hold on;
            count_box=count_box+1;
            %legend(legendInfo,'Location','northwest')
            
            subplot(3,2,5:6)
            p1=plot(ar(rmc>0), rmc(rmc>0),'k:*','DisplayName','odd modes');
            ylabel('Critical 1/\eta')
            xlabel('Box Aspect ratio L_x/L_z')
            hold on;
        else

            subplot(3,2,2)
            plot(rm_entry, rgr_max, 'color',color{count_box},'LineWidth',1.5)
            xlabel('1/\eta');
            ylabel('Real \sigma'); 
            xlim([0 rm_entry(end)])
            ylim([-0.02 0.02]);
            hold on; 
            %legend(legendInfo,'Location','northwest')
            title(['Asymptotic Approximation for a(x) and B(x) for EVEN modes [2,',num2str(n_max),']']);  

            subplot(3,2,4)
            plot(rm_entry, cgr_max, 'color',color{count_box},'LineWidth',1.5)
            xlabel('1/\eta');
            ylabel('Complex \sigma');
            xlim([0 rm_entry(end)])
            hold on;
            count_box=count_box+1; 
            %legend(legendInfo,'Location','northwest')
            
            subplot(3,2,5:6)
            area([0 max(ar)],[0.21 0.21], 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
            text(1, 0.1, 'MEANFIELD', 'color', 'k')
            xlim([0 max(ar)])
            alpha(.02)
            hold on;
            p2=plot(ar(rmc>0), rmc(rmc>0),'b:x','DisplayName','even modes');
            % From Snoopy, critical Rm is far off: 
            ar_snoopy= [10 8];
            rmc_snoopy=[1.1218 1.1458];
            p3=plot(ar_snoopy, rmc_snoopy,'r:o','DisplayName','SNOOPY');
            legend([p1 p2 p3])
           
        end
    end % end m loop 

end
  

