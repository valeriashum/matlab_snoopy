A=1;
B=1;
C=1;
k=1;

files1=[ {'/data/novadisk/vs391/hydro/hydro_lx10_re12/timevar'}...  
        {'/data/novadisk/vs391/hydro/hydro_lx5_re12/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx4_re12/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx2_re12/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx1_re12/timevar'} ];
legendInfo1 = [      {'Box=[10,2,1]'}...
                    {'Box=[5,2,1]'}...
                    {'Box=[4,2,1]'}...
                    {'Box=[2,2,1]'}...
                    {'Box=[1,1,1]'} ];    
                
files2=[{'/data/novadisk/vs391/hydro/hydro_lx10_re12/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re6/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re5/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re4/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re3/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re2/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1p92/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1p5/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1/timevar'}];
legendInfo2 = [     {'R_e=12'}...
                    {'R_e=6'}...
                    {'R_e=5'}...
                    {'R_e=4'}...
                    {'R_e=3'}...
                    {'R_e=2'}...
                    {'R_e=1.92'}...
                    {'R_e=1.5'}...
                    {'R_e=1'} ];  
files3=[{'/data/novadisk/vs391/hydro/hydro_lx1_re12/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx1_re13/timevar'}];
legendInfo3 = [     {'R_e=12'}...
                    {'R_e=13'}];
                
files4=[{'/data/novadisk/vs391/hydro/hydro_lx10_re1/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re2/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re3/timevar'}...
        {'/data/novadisk/vs391/hydro/hydro_lx100_re1/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx100_re2/timevar'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx100_re3/timevar'}];                
legendInfo4 = [     {'[10,2,1] R_e=1'}...
                    {'[10,2,1] R_e=2'}...
                    {'[10,2,1] R_e=3'}...
                    {'[100,2,1] R_e=1'}...
                    {'[100,2,1] R_e=2'}...
                    {'[100,2,1] R_e=3'}];               
files5=[{'/data/novadisk/vs391/hydro/hydro_lx10_re12/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re6/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re5/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re4/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re3/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re2/spectrum.dat'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1p92/spectrum.dat'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1p5/spectrum.dat'}...
        {'/data/novadisk/vs391/hydro/hydro_lx10_re1/spectrum.dat'}];
legendInfo5 = [     {'R_e=12'}...
                    {'R_e=6'}...
                    {'R_e=5'}...
                    {'R_e=4'}...
                    {'R_e=3'}...
                    {'R_e=2'}...
                    {'R_e=1.92'}...
                    {'R_e=1.5'}...
                    {'R_e=1'} ];                 
 files6=[{'/data/novadisk/vs391/hydro/hydro_lx10_re1/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re2/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx10_re3/spectrum.dat'}...
        {'/data/novadisk/vs391/hydro/hydro_lx100_re1/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx100_re2/spectrum.dat'}... 
        {'/data/novadisk/vs391/hydro/hydro_lx100_re3/spectrum.dat'}];                
legendInfo6 = [     {'[10,2,1] R_e=1'}...
                    {'[10,2,1] R_e=2'}...
                    {'[10,2,1] R_e=3'}...
                    {'[100,2,1] R_e=1'}...
                    {'[100,2,1] R_e=2'}...
                    {'[100,2,1] R_e=3'}];                   
sizearray=[size(files1,2) size(files2,2) size(files3,2) size(files4,2) size(files5,2)];    
maxsize=max(sizearray);
for co1=1:maxsize
    color{co1}=rand(1,3); 
end
color{size(files1,2)}=[1,0,0];
for j=1:size(files1,2)
    full_file=importdata(files1{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar);

    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    figure(1)
    if exist('t')&&exist('ev')
        %subplot(3,2,1:2)
        set(gca, 'FontSize', 12)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title('Hydrodynamics of the ABC flow at $R_e=12$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T [T_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K [L^2T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo1,'Location','eastoutside')
%         subplot(3,2,3:4)
%             plot(t,ev-mean(ev),'color',color{j},'LineWidth',1.5);
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             ylabel('$E_K - \bar{E}_K$','fontsize',16, 'Interpreter', 'latex'  );
            
    end

   
%     if exist('t')&&exist('hv')
%         subplot(3,2,5)
%             plot(t,hv,'color',color{j},'LineWidth',1.5);
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             title('$\textit{Helicity } \langle \vec{u}\cdot \vec{\omega} \rangle$','fontsize',16, 'Interpreter', 'latex');
%     end
% 
%     if exist('t')&&exist('w2')
%         subplot(3,2,6)
%             plot(t,w2,'color',color{j},'LineWidth',1.5); 
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             title('$\textit{Vorticity } \int |\vec{\omega}|^2 dV$','fontsize',16, 'Interpreter', 'latex');
%             
%     end
end 
for j=1:size(files2,2)
    full_file=importdata(files2{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar);

    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    figure(2)
    if exist('t')&&exist('ev')
       % subplot(3,2,1:2)
            set(gca, 'FontSize', 12)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title('Hydrodynamics of the ABC flow in a box $[10,2,1]$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T [T_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K [L^2T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo2,'Location','eastoutside')
%         subplot(3,2,3:4)
%             plot(t,ev-mean(ev),'color',color{j},'LineWidth',1.5);
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             ylabel('$E_K - \bar{E}_K$','fontsize',16, 'Interpreter', 'latex'  );
            
    end

   
%     if exist('t')&&exist('hv')
%         subplot(3,2,5)
%             plot(t,hv,'color',color{j},'LineWidth',1.5);
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             title('$\textit{Helicity } \langle \vec{u}\cdot \vec{\omega} \rangle$','fontsize',16, 'Interpreter', 'latex');
%     end
% 
%     if exist('t')&&exist('w2')
%         subplot(3,2,6)
%             plot(t,w2,'color',color{j},'LineWidth',1.5); 
%             hold on;
%             xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
%             title('$\textit{Vorticity } \int |\vec{\omega}|^2 dV$','fontsize',16, 'Interpreter', 'latex');           
%    end
end 
if (1==2)
for j=1:size(files3,2)
    full_file=importdata(files3{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar);

    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    figure(3)
    if exist('t')&&exist('ev')
        subplot(3,2,1:2)
        t(gca, 'FontSize', 12)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title('Hydrodynamics of the ABC flow in a box $[1,1,1]$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo3,'Location','northeast')
        subplot(3,2,3:4)
            plot(t,ev-mean(ev),'color',color{j},'LineWidth',1.5);
            hold on;
            xlabel('$T [T_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K - \bar{E}_K$','fontsize',16, 'Interpreter', 'latex'  );
            
    end

   
    if exist('t')&&exist('hv')
        subplot(3,2,5)
            plot(t,hv,'color',color{j},'LineWidth',1.5);
            hold on;
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('$\textit{Helicity } \langle \vec{u}\cdot \vec{\omega} \rangle$','fontsize',16, 'Interpreter', 'latex');
    end

    if exist('t')&&exist('w2')
        subplot(3,2,6)
            plot(t,w2,'color',color{j},'LineWidth',1.5); 
            hold on;
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('$\textit{Vorticity } \int |\vec{\omega}|^2 dV$','fontsize',16, 'Interpreter', 'latex');
            
    end
end
for j=1:size(files4,2)
    full_file=importdata(files4{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar);

    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    figure(4)
    if exist('t')&&exist('ev')
        subplot(3,2,1:2)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title('$\textit{Hydrodynamics of the ABC flow in a box } [1,1,1]$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo4,'Location','northeast')
        subplot(3,2,3:4)
            plot(t,ev-mean(ev),'color',color{j},'LineWidth',1.5);
            hold on;
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K - \bar{E}_K$','fontsize',16, 'Interpreter', 'latex'  );
            
    end

   
    if exist('t')&&exist('hv')
        subplot(3,2,5)
            plot(t,hv,'color',color{j},'LineWidth',1.5);
            hold on;
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('$\textit{Helicity } \langle \vec{u}\cdot \vec{\omega} \rangle$','fontsize',16, 'Interpreter', 'latex');
    end

    if exist('t')&&exist('w2')
        subplot(3,2,6)
            plot(t,w2,'color',color{j},'LineWidth',1.5); 
            hold on;
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('$\textit{Vorticity } \int |\vec{\omega}|^2 dV$','fontsize',16, 'Interpreter', 'latex');
            
    end
end 
 
for j=1:size(files5,2)
    B_flag = 0;
    u_flag = 1;
    mean_start=1;
    nspec=21; 
    spectruml=importdata(files5{1,j});
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
    
    figure(5)
    if (B_flag==1 && u_flag==1)
        subplot(2,1,1)
            loglog(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
            hold on;
            title('$\textit{Kinetic Energy Spectrum}$','fontsize',20, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
            legend('Kinetic')
             legend(legendInfo4,'Location','northeast')
        subplot(2,1,2)
            loglog(spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
            hold on;
             title('$\textit{Magnetic Energy Spectrum}$','fontsize',20, 'Interpreter', 'latex');
            xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            
    elseif (B_flag==1 && u_flag==0)
        loglog(spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
        hold on;
        title('$\textit{Magnetic Energy Spectrum}$','fontsize',20, 'Interpreter', 'latex');
        xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
        ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
        legend(legendInfo4,'Location','northeast')
    elseif (B_flag==0 && u_flag==1)
        semilogy(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)),...
                    ':o',...
                    'color',color{j},...
                    'LineWidth',1.5,...
                    'MarkerEdgeColor',color{j},...
                    'MarkerFaceColor',color{j},...
                    'MarkerSize',5);
        hold on;
        %title('$\textit{Kinetic Energy Spectrum in a box } [10,2,1]$', 'Interpreter', 'latex','fontsize',20, [0 0 10 10]);
        xlabel('$k \textit{ } [2\pi L^{-1}]$', 'Interpreter', 'latex','fontsize',20);
        ylabel('$E_K \textit{ } [L^2T^{-2}_{turnover}]$', 'Interpreter', 'latex','fontsize',20);
        legend(legendInfo5,'Location','northeast')
   end
end 

for j=1:size(files6,2)
    B_flag = 0;
    u_flag = 1;
    mean_start=1;
    nspec=21; 
    spectruml=importdata(files6{1,j});
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
    
    figure(6)
    if (B_flag==1 && u_flag==1)
        subplot(2,1,1)
            loglog(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)) ...
                    ,'color',color{j},'LineWidth',1.5);
            hold on;
            title('$\textit{Kinetic Energy Spectrum}$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
            legend('Kinetic')
             legend(legendInfo6,'Location','northeast')
        subplot(2,1,2)
            loglog(spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
            hold on;
             title('$\textit{Magnetic Energy Spectrum}$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            
    elseif (B_flag==1 && u_flag==0)
        loglog(spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
        hold on;
        title('$\textit{Magnetic Energy Spectrum}$','fontsize',16, 'Interpreter', 'latex');
        xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
        ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
        legend(legendInfo6,'Location','northeast')
    elseif (B_flag==0 && u_flag==1)
        semilogy(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:))...
                    ,'color',color{j},'LineWidth',1.5);
        hold on;
        title('$\textit{Kinetic Energy Spectrum}$','fontsize',16, 'Interpreter', 'latex');
        xlabel('$k$','fontsize',16, 'Interpreter', 'latex');
        ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
        legend(legendInfo6,'Location','northeast')
    end
end 

end