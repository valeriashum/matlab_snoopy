files1=[ %{'/data/novadisk/vs391/hydro/u_iii/u_iii_no_forcing/timevar'}...  
         %{'/data/novadisk/vs391/hydro/u_iii/u_iii_abc_forcing/timevar'}...  
         {'/data/novadisk/vs391/hydro/u_iii/u_iii_u_iii_forcing/timevar'}...  
         {'/data/novadisk/vs391/hydro/u_iii/u_iii_u_iii_extra_forcing/timevar'}];
legendInfo1 = [ %{'$\vec{f} = 0 $'}...
                %{'$\vec{f} = \nu k^2 \vec{u}_{abc} $'}...
                {'$\vec{f} = \nu k^2 \vec{u}_{iii} $'}...
                {'$\vec{f} = \nu k^2 \vec{u}_{iii} + (\vec{u}_{iii}\cdot \nabla \vec{u}_{iii})^{div-free} $'}];
            
            
for j=1:size(files1,2)
    color{1}=  [1,0,0];
    color{2} = [0,0,0];
    full_file=importdata(files1{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar); 
    
    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    
    figure(1)
    if exist('t')&&exist('ev')
        subplot(3,1,1)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title({'Hydrodynamics of the modified flow at  $R_e=1$ in a box $[10,2,1]$', '\hspace{100pt}Kinetic Energy $E_K$'},'fontsize',18, 'Interpreter', 'latex');
            set(gca, 'FontSize', 12)
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            %ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo1,'Location','east','fontsize',16, 'Interpreter', 'latex') 
            xlim([0 500])
    end
     if exist('t')&&exist('hv')
        subplot(3,1,2)
            plot(t,hv,'color',color{j},'LineWidth',1.5);
            hold on; 
            set(gca, 'FontSize', 12)
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('Helicity $ \langle \vec{u}\cdot \vec{\omega} \rangle$','fontsize',16, 'Interpreter', 'latex');
            xlim([0 500])
    end

    if exist('t')&&exist('w2')
        subplot(3,1,3)
            plot(t,w2,'color',color{j},'LineWidth',1.5); 
            hold on;
            set(gca, 'FontSize', 12)
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            title('Vorticity $\int |\vec{\omega}|^2 dV$','fontsize',16, 'Interpreter', 'latex');
            xlim([0 500])
    end
end
