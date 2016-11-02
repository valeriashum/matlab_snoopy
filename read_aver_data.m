clear all;
rep = ['~/Dropbox/Phd/'];
your_filename = [rep,'by.dat'];

full_file=importdata(your_filename);
nx = full_file(2);
ny = full_file(3);
nz = full_file(4);
gp = full_file(5);

clear full_file your_filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_x = [rep,'bx_yzav.dat'];
filename_y = [rep,'by_yzav.dat'];
filename_z = [rep,'bz_yzav.dat'];

full_file_x=fopen(filename_x,'r');
full_file_y=fopen(filename_y,'r');
full_file_z=fopen(filename_z,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlines = 0;
while ~feof(full_file_x)
    s = fgetl(full_file_x);
    nlines = nlines + 1;
end
numstates = nlines/(nx+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(full_file_x);
fclose(full_file_y);
fclose(full_file_z);
full_file_x=fopen(filename_x,'r');
full_file_y=fopen(filename_y,'r');
full_file_z=fopen(filename_z,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros([nx 1]);
for ind=1:nx
	x(ind)=gp*nx*(ind-1)/(nx-1);
end
clear ind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = zeros([numstates 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bxyzav=zeros(numstates,nx,'like',1 + 1i);
byyzav=zeros(numstates,nx,'like',1 + 1i);
bzyzav=zeros(numstates,nx,'like',1 + 1i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ~feof(full_file_x)
    state = 0;
    for ind=1:numstates
        state = state + 1;
        t(state) = str2num(fgetl(full_file_x));
        t(state) = str2num(fgetl(full_file_y));
        t(state) = str2num(fgetl(full_file_z));
        %fprintf('t=%f\n',t(state,1))
        
        for ind1=1:nx
            ax = str2num(strrep(strrep(strrep(fgetl(full_file_x),'(',''),')',''),',',' '));
            ay = str2num(strrep(strrep(strrep(fgetl(full_file_y),'(',''),')',''),',',' '));
            az = str2num(strrep(strrep(strrep(fgetl(full_file_z),'(',''),')',''),',',' '));
            bxyzav(state,ind1) = ax(1) + ax(2)*1i;
            byyzav(state,ind1) = ay(1) + ay(2)*1i;
            bzyzav(state,ind1) = az(1) + az(2)*1i;
            clear a
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(full_file_x);
fclose(full_file_y);
fclose(full_file_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bxyzav = real(bxyzav)./ max(max(real(bxyzav)));
byyzav = real(byyzav)./ max(max(real(byyzav)));
bzyzav = real(bzyzav)./ max(max(real(bzyzav)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig = figure(1);
set(hFig, 'Position', [100, 50, 1100, 900]); 
set(gca, 'FontSize', 14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1)
    contour(t,x,transpose(real(bxyzav)),...
            'Fill','on');
             xlabel('Time'); ylabel('x'); title('B_x^{norm}');
        colormap jet
        colorbar
        %set(gca,'xscale','log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)
    contour(t,x,transpose(real(byyzav)),...
            'Fill','on');
             xlabel('Time'); ylabel('x'); title('B_y^{norm}');
        colormap jet
        colorbar
        %set(gca,'xscale','log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)
    contour(t,x,transpose(real(bzyzav)),...
            'Fill','on');
             xlabel('Time'); ylabel('x'); title('B_z^{norm}');
        colormap jet
        colorbar
        %set(gca,'xscale','log')
