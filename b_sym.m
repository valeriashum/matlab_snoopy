% Analysing data output from b_symmetry_opt.pro

close all
rep1='/store/ASTRO/vs391/kinematic_dynamo/u_iii';
rep2='/store/ASTRO/vs391/kinematic_dynamo/u_abc';

box=[8 10 16 20 32];

files1=[{[rep1,'/b_sym_8.dat']}...
        {[rep1,'/b_sym_10.dat']}...
        {[rep1,'/b_sym_16.dat']}...
        {[rep1,'/b_sym_20.dat']}...
        {[rep1,'/b_sym_32.dat']} ];   

files2=[{[rep2,'/b_sym_8.dat']}...
        {[rep2,'/b_sym_10.dat']}...
        {[rep2,'/b_sym_16.dat']}...
        {[rep2,'/b_sym_20.dat']}...
        {[rep2,'/b_sym_32.dat']} ];    
    

for j=1:size(files1,2)
    full_file=importdata(files1{1,j});
    timevar=full_file;%.data;
    nvar=size(timevar,2);
    ndat=size(timevar,1);
    % Create a table with all variables
    %tblA = table(Lx,Rm,Gr,lBx,lBy,lBz);
    tblA = table(timevar(:,1),timevar(:,2), timevar(:,3),...
                timevar(:,4), timevar(:,5), timevar(:,6),...
                timevar(:,7), timevar(:,8), timevar(:,9),...
                timevar(:,10),timevar(:,11),timevar(:,12),...
                timevar(:,13),timevar(:,14),timevar(:,15),...
                timevar(:,16),timevar(:,17),timevar(:,18),...
                timevar(:,19),timevar(:,20),timevar(:,21),... 
                timevar(:,22),timevar(:,23),timevar(:,24),... 
                timevar(:,25));
    % Sort the rows of the table based on Rm
    tblB = sortrows(tblA,1); 
    Rm      = tblB{:,1}; 
    % x0 is the shift off x=0 for optimal symmetry in x-direction   
    x0_x_xy = tblB{:,2}; 
    x0_y_xy = tblB{:,3}; 
    x0_z_xy = tblB{:,4}; 
    x0_x_xz = tblB{:,5}; 
    x0_y_xz = tblB{:,6}; 
    x0_z_xz = tblB{:,7}; 
    
    % best symmetry in x-direction for two halfs
    d0_x_xy = tblB{:,8}; 
    d0_y_xy = tblB{:,9}; 
    d0_z_xy = tblB{:,10}; 
    d0_x_xz = tblB{:,11}; 
    d0_y_xz = tblB{:,12}; 
    d0_z_xz = tblB{:,13}; 
    % mean(abs(B^I - B^II)) for the symmetrized case 
    d_x_xy  = tblB{:,14}; 
    d_y_xy  = tblB{:,15}; 
    d_z_xy  = tblB{:,16}; 
    d_x_xz  = tblB{:,17}; 
    d_y_xz  = tblB{:,18}; 
    d_z_xz  = tblB{:,19}; 
    % shift need for symmatrisation
    s_x_xy  = tblB{:,20}; 
    s_y_xy  = tblB{:,21}; 
    s_z_xy  = tblB{:,22}; 
    s_x_xz  = tblB{:,23}; 
    s_y_xz  = tblB{:,24}; 
    s_z_xz  = tblB{:,25};  
    
      
hFig = figure(j);
set(hFig, 'Position', [100*2*j-100, 50, 400, 400]);
    % XXXXXXXX
    subplot(4,3,1)
        plot(Rm, x0_x_xy/(6.28),'ro',Rm, x0_x_xz/(6.28),'k*') 
            title('B_x')
            ylabel('x_0')

        ylim([0 box(j)])
    subplot(4,3,4)
        plot(Rm,d0_x_xy/(6.28) ,'ro',Rm,d0_x_xz/(6.28) ,'k*') 
            ylabel('S_x')
 
    subplot(4,3,7)
        plot(Rm,s_x_xy/(6.28) ,'ro',Rm,s_x_xz/(6.28) ,'k*') 
            ylabel('y_0, z_0')

        ylim([0 4])
    subplot(4,3,10)
        plot(Rm,d_x_xy, 'ro',Rm,d_x_xz ,'k*') 
            ylabel('S_y,S_z')
            xlabel('L_x')
        ylim([0 0.4])
   % YYYYYY
    subplot(4,3,2)
        plot(Rm, x0_y_xy/(6.28),'ro',Rm, x0_y_xz/(6.28),'k*') 
            title('B_y')
        ylim([0 box(j)])
    subplot(4,3,5)
        plot(Rm,d0_y_xy/(6.28) ,'ro',Rm,d0_y_xz/(6.28) ,'k*') 
        
    subplot(4,3,8)
        plot(Rm,s_y_xy/(6.28) ,'ro',Rm,s_y_xz/(6.28) ,'k*') 
        ylim([0 4])
    subplot(4,3,11)
        plot(Rm,d_y_xy, 'ro',Rm,d_y_xz ,'k*') 
            xlabel('L_x')
             ylim([0 0.4])
      % ZZZZZZ
    subplot(4,3,3)
        plot(Rm, x0_z_xy/(6.28),'ro',Rm, x0_z_xz/(6.28),'k*') 
            title('B_z')
         ylim([0 box(j)])
    subplot(4,3,6)
        plot(Rm,d0_z_xy/(6.28) ,'ro',Rm,d0_z_xz/(6.28) ,'k*') 
        
    subplot(4,3,9)
        plot(Rm,s_z_xy/(6.28) ,'ro',Rm,s_z_xz/(6.28) ,'k*') 
        ylim([0 4])
    subplot(4,3,12)
        plot(Rm,d_z_xy, 'ro',Rm,d_z_xz ,'k*') 
              xlabel('L_x')
          ylim([0 0.4])
        
        
        
end