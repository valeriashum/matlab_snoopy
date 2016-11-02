close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep=[{'/store/ASTRO/vs391/nonlinear_dynamo/u_abc'     , '/store/ASTRO/vs391/nonlinear_dynamo/u_iii'}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box=[10 16 32];
for ii=1:size(rep,2)
    clear tblA tblB run cas 
    for jj=1:size(box,2)
        rep1=[rep{ii}, num2str(box(jj))];
        file_info=[rep1,'/info'];
        full_file=importdata(file_info);
        info=full_file;
        tblA = table(info(:,1),info(:,2), info(:,3), info(:,4));
        % Sort the rows of the table based on Re, then Rm
        tblB = sortrows(tblA,[4,3]); 
        run = tblB{1:end,1}; 
        cas = tblB{1:end,2}; 
        etainv = tblB{1:end,3};  
        nuinv  = tblB{1:end,4}; 
        for j=1:size(nuinv,1)
            clear x t bxyzav byyzav bzyzav
            if run(j) < 100
                rep2 = [rep,'/run_',num2str(run(j)),'/for_analysis/graphics/'];
            else
                rep2 = [rep,'/run_',num2str(run(j)),'/for_analysis/result1/'];
            end
            your_filename = [rep2,'by.dat'];

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
            hFig = figure(10*ii + jj);
            set(hFig, 'Position', [100, 50, 1100, 900]); 
            set(gca, 'FontSize', 14)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rows = 4;
            if mod(size(Ree,1),rows)==0
                cols = fix(size(etainv,1)/rows);
            else 
                cols = fix(size(etainv,1)/rows) + 1;
            end
            subplot(cols,rows,j)
                contour(t,x,transpose(real(bxyzav)),...
                        'Fill','on');
                         xlabel('Time'); ylabel('x'); title('B_x^{norm}');
                    colormap jet
                    colorbar
        end
    end
end
    
        