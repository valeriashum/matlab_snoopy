
files2=[{'/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/results/box_8/from_summary'}... 
        {'/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/results/box_10/from_summary'}...
        {'/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/results/box_16/from_summary'}];
titleInfo2 = [     {'[8,2,1]'}...
                    {'[10,2,1]'}...
                    {'[16,2,1]'}];  
          
 color{3} = [0.33,0.33,1.0]; color{4} = [0,0.5,0]; color{5} = [0.5,0,0.5];                
 for j=1:size(files2,2)
        full_file=importdata(files2{1,j});
        timevar=full_file.data;
        nvar=size(timevar,2);
        n=size(timevar,1);
        % Create a table with all variables
        % Box Rm  Gr  KX  KY  KZ  EM  KX1 KY1 KZ1 EM1 LBX LBY LBZ
        %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
        tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6), ...
                     timevar(:,7),timevar(:,8),timevar(:,9),timevar(:,10),timevar(:,11),timevar(:,12), ...
                     timevar(:,13),timevar(:,14));
        clear timevar;
        % Sort the rows of the table based on Rm
        tblB = sortrows(tblA,2);
        Lx = tblB{:,1};
        Rm = tblB{:,2};
        Gr = tblB{:,3};
        KX = tblB{:,4};
        KY = tblB{:,5}; 
        KZ = tblB{:,6};
        EM = tblB{:,7};
        KX1 = tblB{:,8};
        KY1 = tblB{:,9};
        KZ1 = tblB{:,10};
        EM1 = tblB{:,11};
        lBx= tblB{:,12};
        lBy= tblB{:,13};
        lBz= tblB{:,14};
        
        legendInfo = [{'dominant k_x'} {'sub-dominant k_x'} {'(l_B)^{-1}'}];      
        maekrsize=15.0;
        
        subplot(1,size(files2,2),j)
        plot(Rm,KX, 'o', 'color',color{j+2},'MarkerSize',maekrsize );
        title(titleInfo2(j));
        xlabel('R_m');
        ylabel('k_x');
         set(gca, 'FontSize', 14)     
        hold on;
        plot(Rm, KX1,'x', 'color',color{j+2},'MarkerSize',maekrsize);
        plot(Rm, 1./(lBy/6.28),'*', 'color',color{j+2},'MarkerSize',maekrsize);
        xlim([0,2])
        ylim([0.05,0.45])
        legend(legendInfo,'Location','southoutside','fontsize',16, 'Box', 'off');
    
            for ind=1./Lx(1):1./Lx(1):0.45
                plot(xlim,[ind ind],':', 'color',color{j+2} );
            end
       
 end 