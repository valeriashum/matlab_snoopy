clear all;
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
for i=1:1000
    ar(i) = i*5;    
end

for co=1:2
    n_max=nn(co);
    count_box=1;
    rmc = zeros(1,size(ar,2));
    for box = 1:size(ar,2) 
        color{box} = rand(1,3);
        % wavenumber corresponding to the longest wavelength on the box
        m =1./ar(box);

        count=1; % data output count for each rm
        % magnetic diffusivity & magnetic Reynolds number
        for rm=0.0001:0.0001:2.0
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
            %fprintf('box=%d,rm=%f,m=%f,ar(box), D(1,1)=%f\n',box, rm,m,ar(box),D(1,1)); 
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
             
        % Critical Rmc 
            [min_gr, II] = min(rgr_max(rgr_max>0));
            rm_entry_positive = rm_entry(rgr_max>0);
            rmc(box) = rm_entry_positive(II);
            if (mod(n_max,2)~=0)
                ar_odd = ar(rmc>0);
                rmc_odd= rmc(rmc>0);
            else
                ar_even = ar(rmc>0);
                rmc_even= rmc(rmc>0);
            end
      end % end m loop 
end
%% 
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
  
figure(3)
    semilogx(ar, rmc_odd.*sqrt(ar),':r', ar, rmc_even.*sqrt(ar),':b', plot(ar(rmc_odd>0), rmc_even(rmc_odd>0)/rmc_odd(rmc_odd>0),'ko')) 
    xlabel('L_x [L]');
    legend('R_{mc}^{odd}/R_{mc}^{abc}', 'R_{mc}^{even}/R_{mc}^{abc}','R_{mc}^{even}/R_{mc}^{odd}','fontsize',18, 'Interpreter', 'latex' );
        