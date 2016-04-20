% Get the data
% @i is only used for outputting to the table
% save this file in KinematicOutput/graphics
rep ='/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/kinematicOutput_box_10_3/';
rep1='/media/vs391/TOSHIBAHARD/WORK/u_abc/box_10_kinematicOutput_5/';

for constt = 1  % file number
         fprintf('Begin: modes%d\n', constt);
         clearvars -except const constt folder rep rep1
         
         filename_in = [rep, '/rates' num2str(constt) '/modes' num2str(constt)];
         filename_gr = [rep, 'timevar' num2str(constt) '.dat'];
        try
            % get the magnetic Reynolds number for the sim
            file=importdata(filename_in);       file_gr = importdata(filename_gr);
            timevar=file.data;                  timevar_gr = file_gr.data;
            nvar=size(timevar,2);               nvar_gr=size(timevar_gr,2);

            % Create a table with all variables
            %Box	Rm	KX	KY	KZ	RGR	RSq_r	Bntfinnorm	CGR	Period
            tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),...
                         timevar(:,6),timevar(:,7),timevar(:,8),timevar(:,9),timevar(:,10));
            % Sort the rows of the table
            tblB = sortrows(tblA,6);
            
            % Save the variables (10th entry is arbitrary)
            gr  = timevar_gr(end,nvar_gr);
            box = timevar(10,1);
            Rm  = timevar(10,2);
            KX  = timevar(:,3);
            KY  = timevar(:,4);
            KZ  = timevar(:,5);
            RGR = timevar(:,6);
            Bn  = timevar(:,8);
            
            % Calculate mean and median of the velocity field
            meany   = mean(sort(RGR(RGR < gr & RGR > 0)));
            mediany = median(sort(RGR(RGR < gr & RGR > 0))); 
            
            meany_long = mean(RGR(KX <= 0.5));
            mediany_long = median(RGR(KX <= 0.5));
            
            meany_short = mean(RGR(KX > 0.5));
            mediany_short = median(RGR(KX > 0.5));
            
            % Calculate power at each kx = mod of 1, 0.1, 0.9
            gr_kx1_ky0 = zeros(size(KX,1),1);           gr_kx0p9_ky0 = zeros(size(KX,1),1);
            gr_kx1_ky0p5 = zeros(size(KX,1),1);         gr_kx0p9_ky0p5= zeros(size(KX,1),1);
            gr_kx0p1_ky0 = zeros(size(KX,1),1);         gr_kx1p1_ky0= zeros(size(KX,1),1);
            gr_kx0p1_ky0p5 = zeros(size(KX,1),1);       gr_kx1p1_ky0p5= zeros(size(KX,1),1);
            
            for ii=1:size(KX,1)
                if (RGR(ii) < gr) && (RGR(ii) > min(meany,mediany));
                    if mod(abs(KX(ii))-round(abs(KX(ii))),1) < 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),1)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx1_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),0.1)< 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx0p1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),0.1)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx0p1_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),0.9) < 1.e-4&& mod(abs(KY(ii)),1)==0
                        gr_kx0p9_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),0.9)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx0p9_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),1.1) < 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx1p1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),1.1) < 1.e-4 && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx1p1_ky0p5(ii,1) = RGR(ii);
                    end
                end
            end
            legendInfo= [ {'k_f=(0.1,0,0)'}...
                         {'k_f=(0.1,0.5,0)'}...
                         {'k_f=(1,0,0)'} ...
                         {'k_f=(1,0.5,0)'} ...
                         {'k_f=(1.1,0,0)'}...
                         {'k_f=(1.1,0.5,0)'}...
                         {'k_f=(0.9,0,0)'}...
                         {'k_f=(0.9,0.5,0)'}  ];
            
           
            figure(1)
            subplot(2,2,1) 
                plot(sort(RGR(RGR < gr & RGR > 0)), '-k', 'Linewidth', 1.5);
                hold on;
                plot(get(gca,'xlim'), [meany meany] ,':r', 'Linewidth', 2.5);
                plot(get(gca,'xlim'), [mediany mediany],':b', 'Linewidth', 2.5);
                plot(get(gca,'xlim'), [meany_long meany_long] ,'--r', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [mediany_long mediany_long],'--b', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [meany_short meany_short] ,'-.r', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [mediany_short mediany_short],'-.b', 'Linewidth', 1.5);
                xlabel('n^{th} mode');
                ylabel('\sigma (k)');  
                legend([{'\sigma'} {'mean'} {'median'} {'mean_{k_x < 0.5}'} {'median_{k_x < 0.5}'} {'mean_{k_x > 0.5}'} {'median_{k_x > 0.5}'}],'Location','southeast');
                title([{ '\hspace{400pt} Kinematic dynamo in the presence of $u_{iii}$ at $R_m = 0.7$'}, {'\\'}],'fontsize',16, 'Interpreter', 'latex');
                
            subplot(2,2,2)
                plot(sqrt(KX.^2+KY.^2+KZ.^2),Bn/sqrt(sum(Bn.^2)),'ko');
                hold on; 
                plot(sqrt(KX(RGR < gr & RGR > min(meany,mediany)).^2+KY(RGR < gr & RGR > min(meany,mediany)).^2+KZ(RGR < gr & RGR > min(meany,mediany)).^2),Bn(RGR < gr & RGR > min(meany,mediany))/sqrt(sum(Bn(RGR < gr & RGR > min(meany,mediany)).^2)),'ro'); 
                xlabel('|K|');
                ylabel('Power B_n(K)/B_{rms}');   
                legend([{'\sigma < \sigma_{mm}'} {' \sigma_{mm}< \sigma < \sigma_{max}'} ],'Location','northeast', 'fontsize', 14);
                
            subplot(2,2,3)
                plot(sqrt(KX(RGR < gr & RGR > min(meany,mediany) ).^2+KY(RGR < gr & RGR > min(meany,mediany)).^2+KZ(RGR < gr & RGR > min(meany,mediany)).^2),RGR(RGR < gr & RGR >min(meany,mediany)  ),...
                    'ko'); 
                xlabel('|K|');
                ylabel(' \sigma_{mm}< \sigma < \sigma_{max}');   
            subplot(2,2,4)
                
                plot(abs(KX(gr_kx0p1_ky0 > 0 & gr_kx0p1_ky0 < gr)),gr_kx0p1_ky0(gr_kx0p1_ky0 > 0 & gr_kx0p1_ky0 < gr) ,'ro'); 
                hold on;
                plot(abs(KX(gr_kx0p1_ky0p5 > 0 & gr_kx0p1_ky0p5 < gr)),gr_kx0p1_ky0p5(gr_kx0p1_ky0p5 > 0 & gr_kx0p1_ky0p5 < gr) ,'rx');
                plot(abs(KX(gr_kx1_ky0 > 0 & gr_kx1_ky0 < gr)),gr_kx1_ky0(gr_kx1_ky0 > 0 & gr_kx1_ky0 < gr) ,'ko'); 
                plot(abs(KX(gr_kx1_ky0p5 > 0 & gr_kx1_ky0p5 < gr)),gr_kx1_ky0p5(gr_kx1_ky0p5 > 0 & gr_kx1_ky0p5 < gr) ,'kx'); 
                plot(abs(KX(gr_kx0p9_ky0 > 0 & gr_kx0p9_ky0 < gr)),gr_kx0p9_ky0(gr_kx0p9_ky0 > 0 & gr_kx0p9_ky0 < gr) ,'bo'); 
                plot(abs(KX(gr_kx0p9_ky0p5 > 0 & gr_kx0p9_ky0p5 < gr)),gr_kx0p9_ky0p5(gr_kx0p9_ky0p5 > 0 & gr_kx0p9_ky0p5 < gr) ,'bx');
                plot(abs(KX(gr_kx1p1_ky0 > 0 & gr_kx1p1_ky0 < gr)),gr_kx1p1_ky0(gr_kx1p1_ky0 > 0 & gr_kx1p1_ky0 < gr) ,'go'); 
                plot(abs(KX(gr_kx1p1_ky0p5 > 0 & gr_kx1p1_ky0p5 < gr)),gr_kx1p1_ky0p5(gr_kx1p1_ky0p5 > 0 & gr_kx1p1_ky0p5 < gr) ,'gx');
                xlabel('k_x');
                ylabel(' \sigma_{mm}< \sigma < \sigma_{max}');       
                legend(legendInfo,'Location','northeast');

                        
            catch ME
            fprintf('Error with the modes %d file %s\n',constt, ME.message);
            continue; 
        end   
end
for constt = 2  % file number
         fprintf('Begin: modes%d\n', constt);
         clearvars -except const constt folder rep rep1
         
         filename_in = [rep1, '/rates' num2str(constt) '/modes' num2str(constt)];
         filename_gr = [rep1, 'timevar' num2str(constt) '.dat'];
        try
            % get the magnetic Reynolds number for the sim
            file=importdata(filename_in);       file_gr = importdata(filename_gr);
            timevar=file.data;                  timevar_gr = file_gr.data;
            nvar=size(timevar,2);               nvar_gr=size(timevar_gr,2);

            % Create a table with all variables
            %Box	Rm	KX	KY	KZ	RGR	RSq_r	Bntfinnorm	CGR	Period
            tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),...
                         timevar(:,6),timevar(:,7),timevar(:,8),timevar(:,9),timevar(:,10));
            % Sort the rows of the table
            tblB = sortrows(tblA,6);
            
            % Save the variables (10th entry is arbitrary)
            gr  = timevar_gr(end,nvar_gr);
            box = timevar(10,1);
            Rm  = timevar(10,2);
            KX  = timevar(:,3);
            KY  = timevar(:,4);
            KZ  = timevar(:,5);
            RGR = timevar(:,6);
            Bn  = timevar(:,8);
            plot(KX, RGR)
            % Calculate mean and median of the velocity field
            meany   = mean(sort(RGR(RGR < gr & RGR > 0)));
            mediany = median(sort(RGR(RGR < gr & RGR > 0))); 
            
            meany_long = mean(RGR(KX <= 0.5));
            mediany_long = median(RGR(KX <= 0.5));
            
            meany_short = mean(RGR(KX > 0.5));
            mediany_short = median(RGR(KX > 0.5));
            
            % Calculate power at each kx = mod of 1, 0.1, 0.9
            gr_kx1_ky0 = zeros(size(KX,1),1);           gr_kx0p9_ky0 = zeros(size(KX,1),1);
            gr_kx1_ky0p5 = zeros(size(KX,1),1);         gr_kx0p9_ky0p5= zeros(size(KX,1),1);
            gr_kx0p1_ky0 = zeros(size(KX,1),1);         gr_kx1p1_ky0= zeros(size(KX,1),1);
            gr_kx0p1_ky0p5 = zeros(size(KX,1),1);       gr_kx1p1_ky0p5= zeros(size(KX,1),1);
            
            for ii=1:size(KX,1)
                if (RGR(ii) < gr) && (RGR(ii) > min(meany,mediany));
                    if mod(abs(KX(ii))-round(abs(KX(ii))),1) < 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),1)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx1_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),0.1)< 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx0p1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii))-round(abs(KX(ii))),0.1)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx0p1_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),0.9) < 1.e-4&& mod(abs(KY(ii)),1)==0
                        gr_kx0p9_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),0.9)< 1.e-4  && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx0p9_ky0p5(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),1.1) < 1.e-4 && mod(abs(KY(ii)),1)==0
                        gr_kx1p1_ky0(ii,1) = RGR(ii);
                    end
                    if mod(abs(KX(ii)),1.1) < 1.e-4 && mod(abs(KY(ii))-round(abs(KY(ii))),0.5)< 1.e-4
                        gr_kx1p1_ky0p5(ii,1) = RGR(ii);
                    end
                end
            end
            legendInfo= [ {'k_f=(1,0,0)'} ...
                         {'k_f=(1,0.5,0)'}  ];
            
           
            figure(2)
            subplot(2,2,1) 
                plot(sort(RGR(RGR < gr & RGR > 0)), '-k', 'Linewidth', 1.5);
                hold on;
                plot(get(gca,'xlim'), [meany meany] ,':r', 'Linewidth', 2.5);
                plot(get(gca,'xlim'), [mediany mediany],':b', 'Linewidth', 2.5);
                plot(get(gca,'xlim'), [meany_long meany_long] ,'--r', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [mediany_long mediany_long],'--b', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [meany_short meany_short] ,'-.r', 'Linewidth', 1.5);
                plot(get(gca,'xlim'), [mediany_short mediany_short],'-.b', 'Linewidth', 1.5);
                xlabel('n^{th} mode');
                ylabel('\sigma (k)');  
                legend([{'\sigma'} {'mean'} {'median'} {'mean_{k_x < 0.5}'} {'median_{k_x < 0.5}'} {'mean_{k_x > 0.5}'} {'median_{k_x > 0.5}'}],'Location','southeast');
                title([{ '\hspace{400pt} Kinematic dynamo in the presence of $u_{abc}$ at $R_m = 0.7$'}, {'\\'}],'fontsize',16, 'Interpreter', 'latex');
                
            subplot(2,2,2)
                plot(sqrt(KX.^2+KY.^2+KZ.^2),Bn/sqrt(sum(Bn.^2)),'ko');
                hold on; 
                plot(sqrt(KX(RGR < gr & RGR > min(meany,mediany)).^2+KY(RGR < gr & RGR > min(meany,mediany)).^2+KZ(RGR < gr & RGR > min(meany,mediany)).^2),Bn(RGR < gr & RGR > min(meany,mediany))/sqrt(sum(Bn(RGR < gr & RGR > min(meany,mediany)).^2)),'ro'); 
                xlabel('|K|');
                ylabel('Power B_n(K)/B_{rms}');   
                legend([{'\sigma < \sigma_{mm}'} {' \sigma_{mm}< \sigma < \sigma_{max}'} ],'Location','northeast', 'fontsize', 14);
                
            subplot(2,2,3)
                plot(sqrt(KX(RGR < gr & RGR > min(meany,mediany) ).^2+KY(RGR < gr & RGR > min(meany,mediany)).^2+KZ(RGR < gr & RGR > min(meany,mediany)).^2),RGR(RGR < gr & RGR >min(meany,mediany)  ),...
                    'ko'); 
                xlabel('|K|');
                ylabel(' \sigma_{mm}< \sigma < \sigma_{max}');   
            subplot(2,2,4)            
                plot(abs(KX(gr_kx1_ky0 > 0 & gr_kx1_ky0 < gr)),gr_kx1_ky0(gr_kx1_ky0 > 0 & gr_kx1_ky0 < gr) ,'ko'); 
                hold on;
                plot(abs(KX(gr_kx1_ky0p5 > 0 & gr_kx1_ky0p5 < gr)),gr_kx1_ky0p5(gr_kx1_ky0p5 > 0 & gr_kx1_ky0p5 < gr) ,'kx'); 
       
                xlabel('k_x');
                ylabel(' \sigma_{mm}< \sigma < \sigma_{max}');       
                legend(legendInfo,'Location','northeast');

                        
            catch ME
            fprintf('Error with the modes %d file %s\n',constt, ME.message);
            continue; 
        end   
end
