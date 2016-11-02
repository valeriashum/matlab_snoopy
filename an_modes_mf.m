close all
rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

files1=[{[rep1,'/modes_10.dat']}... 
        {[rep1,'/modes_10.dat']}...
        {[rep1,'/modes_16.dat']}...  
        {[rep1,'/modes_20.dat']}... 
        {[rep1,'/modes_32.dat']}];
        

files2=[{[rep2,'/modes_8.dat']}... 
        {[rep2,'/modes_10.dat']}...
        {[rep2,'/modes_16.dat']}...  
        {[rep2,'/modes_20.dat']}... 
        {[rep2,'/modes_32.dat']}];
    
    
legendInfo1 = [     {'Box=[8,2,1]'}...
                    {'Box=[10,2,1]'}...
                    {'Box=[16,2,1]'}...
                    {'Box=[20,2,1]'}...
                    {'Box=[32,2,1]'}];  
 s = {'+','o','*','.','x','s','d','^','v','>','<','p','h','+','o','*','.','x','s','d','^','v','>','<','p','h'};  
 for ind=1:size(files1,2)
    color{ind}=[0,0,0];
    color1{ind}=[0.83,0,0];
 end      

for box=1:size(files1,2)

    full_file=importdata(files1{1,box});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6),timevar(:,7));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,[4 3 5 2]);
    
    % Find how many different modes are present
    % @mode_box - array that contains the location of first mode of the group
    % @n_modes  - number of various dynamo modes
    iii=0;
    k=0;
    mode_box{1}=1;
                    % count number of single mode data points at the end
                    % to calculate the total number of modes 
                    % @n_sin is number of single modes
                    n_sin = 0;

                    for j=1:ndat-1
                        if (tblB{j,3}==tblB{j+1,3} & tblB{j,4}==tblB{j+1,4} & tblB{j,5}==tblB{j+1,5})
                            iii=iii+1;
                            n_sin = 0;
                        else
                            n_sin = n_sin+1;    
                            k = k+1;
                            mode_box{k+1}=j+1;
                        end
                    end

                    % For purpose of using mode_box{n_modes + 1} in calculations
                    % Set an extra entry for mode_box = last mode data point
                    mode_box{k+2}=ndat;

                    % Calculate total number of modes growing in the system
                    % Should be equal to size(mode_box,2)-1
                    n_modes = ndat - iii - n_sin;
                    county=1;
                    county1=1;
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%    LOOP    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for i=1:n_modes
                        
                       clear Lx Rm kx ky kz Gr;                      
                       % Get the EM and time data to fit an exp
                       % @x is time
                       % @y is B(k)
                       for j=1:(mode_box{i+1}-mode_box{i})
                            Lx(j) = tblB{j+mode_box{i}-1,1};
                            Rm(j) = tblB{j+mode_box{i}-1,2}.*(1.5+1./Lx(1).^2)/sqrt(1.5*(3+2*1./Lx(1).^4+9*1./Lx(1).^2));
                            kx(j) = tblB{j+mode_box{i}-1,3};
                            ky(j) = tblB{j+mode_box{i}-1,4};
                            kz(j) = tblB{j+mode_box{i}-1,5}; 
                            Gr(j) = tblB{j+mode_box{i}-1,6};
                            Bn(j) = tblB{j+mode_box{i}-1,7};
                       end                       
                       if (kx(1)>0 && ky(1)==0)
                            legendinfo{county}= [num2str(kx(1)) ',' num2str(ky(1)) ',' num2str(kz(1))];
                            county=county+1;
                       end
                       if (kx(1)>0 && ky(1)>0)
                            legendinfo1{county1}= [num2str(kx(1)) ',' num2str(ky(1)) ',' num2str(kz(1))];
                            county1=county1+1;
                       end
                       
                       hFig=figure(box);
                       set(hFig, 'Position',[0 100 1000 600])
                       %subplot(1,size(files1,2),box)
                       subplot(1,3,1)
                            set(gca, 'FontSize', 18)    
                            plot(Rm(Gr>0 & kx>0 & ky==0), Gr(Gr>0 & kx>0 & ky==0),s{county},...
                                    'LineStyle', ':',...
                                    'LineWidth',1.5,...
                                    'MarkerSize',6.5,...
                                    'color',color{box},...
                                    'MarkerEdgeColor', color{box}, ...
                                    'MarkerFaceColor', color{box});
                            hold on;    
                            plot(Rm(Gr>0 & kx>0 & ky>0), Gr(Gr>0 & kx>0 & ky>0),s{county1},...
                                    'LineStyle', ':',...
                                    'LineWidth',1.5,...
                                    'MarkerSize',4.5,...
                                    'color',color1{box},...
                                    'MarkerEdgeColor', color1{box}, ...
                                    'MarkerFaceColor', color1{box}); 
                            xlim([0 1])   
                            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                            ylabel('$\sigma [T^{-1}_{turnover}]$','fontsize',16, 'Interpreter', 'latex');
                        subplot(1,3,2:3)        
                            set(gca, 'FontSize', 18)    
                            plot(Rm(Gr>0 & kx>0 & ky==0), Gr(Gr>0 & kx>0 & ky==0),s{county},...
                                    'LineStyle', ':',...
                                    'LineWidth',1.5,...
                                    'MarkerSize',6.5,...
                                    'color',color{box},...
                                    'MarkerEdgeColor', color{box}, ...
                                    'MarkerFaceColor', color{box});
                            hold on;    
                            plot(Rm(Gr>0 & kx>0 & ky>0), Gr(Gr>0 & kx>0 & ky>0),s{county1},...
                                    'LineStyle', ':',...
                                    'LineWidth',1.5,...
                                    'MarkerSize',4.5,...
                                    'color',color1{box},...
                                    'MarkerEdgeColor', color1{box}, ...
                                    'MarkerFaceColor', color1{box});    
                            title(legendInfo1{box},'fontsize',16, 'Interpreter', 'latex')        
                            xlabel('$R_m$','fontsize',16, 'Interpreter', 'latex');
                            
                            
                    end
C = [legendinfo, legendinfo1];
legend(C, 'Location', 'eastoutside', 'FontSize', 14)
                                     
end





