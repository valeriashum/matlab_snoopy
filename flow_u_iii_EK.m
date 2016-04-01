EKabc= [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 ];
urmsabc=sqrt(2.*EKabc);
 
files=['/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/u_iii_ek_fn_m'];
full_file=importdata(files);


timevar=transpose(full_file.data);

nvar=size(timevar,1);

Lx= timevar(1,:);   
Ly= timevar(2,:);   
Lz= timevar(3,:);   
EK= timevar(4,:);   

m=Lz./Lx;   
urms  = sqrt(2.*EK);
EK = 2*EK;
EKabc = 2*EKabc;
color{1} = [0,0,0];
color{2} = [1,0,0];
figure(1) 
subplot(2,2,1)
    plot(m,EK,':o','LineWidth',1.5, ...
            'color', color{2}, ...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',3.5);
    hold on;
    title(' Kinetic Energy ','fontsize',18, 'Interpreter', 'latex');
   
    ylabel('$<u_iu_i> [L^2T^{-2}_{turnover}]$','fontsize',18, 'Interpreter', 'latex');
    set(gca, 'FontSize', 10) 
    legendInfo1 = [{'$u_{III}$'}...
                   {'$~ m^2+\frac{3}{2}$ '}];
    legend(legendInfo1,'Location','northwest','fontsize',14, 'Interpreter', 'latex'); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% POLYFIT %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % p(x)=p1xn+p2xn−1+...+pnx+pn+1.
    % @f(1) the slope on line fit for n=1
    f1 = polyfit(m,EK,2);
    plot(m,f1(1).*m.^2+f1(2).*m + f1(3), 'b:','LineWidth',1.0)
subplot(2,2,3)    
    plot(m,EKabc, ':o','LineWidth',1.5, ...
            'color', color{1}, ...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',3.5);
legendInfo2 = [{'$u_{ABC} $ '}];         
legend(legendInfo2,'Location','northwest','fontsize',14, 'Interpreter', 'latex');     
xlabel('$m [2\pi L^{-1}]$','fontsize',18, 'Interpreter', 'latex');
set(gca, 'FontSize', 10) 
files1=[ {'/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/snoopy1/kinematicOutput_box_10_1/spectrum1.dat'}];
    
    
legendInfo1 = [     {'u_{abc}'}...
                    {'u_{III}'}];
    mean_start=1;
    nspec=21;    % 9 for old spectrum (no transfer). 15 for recent versions, 21 for v5.0 version
    spectruml=importdata(files1{1,1});
    spectrum.k=spectruml(1,1:(end-1));
    spectrum.n=spectruml(2,1:(end-1));
    spectrum.vx=spectruml(3:nspec:end,2:end);
    spectrum.vy=spectruml(4:nspec:end,2:end);
    spectrum.vz=spectruml(5:nspec:end,2:end);
    spectrum.bx=spectruml(6:nspec:end,2:end);
    spectrum.by=spectruml(7:nspec:end,2:end);
    spectrum.bz=spectruml(8:nspec:end,2:end);
    spectrum.th=spectruml(9:nspec:end,2:end);
    spectrum.vxvy=spectruml(10:nspec:end,2:end);
    spectrum.bxby=spectruml(11:nspec:end,2:end);
    spectrum.ad_vx=spectruml(12:nspec:end,2:end);
    spectrum.ad_vy=spectruml(13:nspec:end,2:end);
    spectrum.ad_vz=spectruml(14:nspec:end,2:end);
    spectrum.tr_bx=spectruml(15:nspec:end,2:end);
    spectrum.tr_by=spectruml(16:nspec:end,2:end);
    spectrum.tr_bz=spectruml(17:nspec:end,2:end);
    spectrum.tr_vx=spectruml(18:nspec:end,2:end);
    spectrum.tr_vy=spectruml(19:nspec:end,2:end);
    spectrum.tr_vz=spectruml(20:nspec:end,2:end);
    spectrum.hel=spectruml(21:nspec:end,2:end)+spectruml(22:nspec:end,2:end)+spectruml(23:nspec:end,2:end);


    spectrum.t=transpose(spectruml(4:nspec:end,1));

    clear spectrum1;
    
    spk =[1.0];
    meanek=[1.5];
    
    subplot(2,2,[2,4]) 
        ax=gca;
        ax.YAxisLocation = 'right';
        semilogy(spectrum.k, mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:))...
            ,':o','LineWidth',1.5, ...
            'color', color{2}, ...
            'MarkerEdgeColor', color{2}, ...
            'MarkerFaceColor', color{2}, ...
            'MarkerSize',3.5);
        ylim([1.e-30 2.0]);
        title('Kinetic Energy Spectrum','fontsize',18, 'Interpreter', 'latex');
        xlabel('$k [2\pi L^{-1}]$','fontsize',16, 'Interpreter', 'latex');
        ylabel('$E_K [L^2T_{turnover}^{-2}]$','fontsize',16, 'Interpreter', 'latex'  );
        set(gca, 'FontSize', 10) 
         
         %legend(legendInfo1,'Location','northeast')
        
        hold on;
        semilogy(spk,meanek ...
            ,'o','LineWidth',1.5, ...
            'color', color{1}, ...
            'MarkerEdgeColor', color{1}, ...
            'MarkerFaceColor', color{1}, ...
            'MarkerSize',3.5);





