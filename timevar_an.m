
% Example File: how to read the timevar file in Matlab
% This example reads all the variable in the timevar file
% and display the time history of several quantities...

% USE @constt TO LOOP OVER
% MULTIPLE TIMEVAR FILES
for constt =0:15
    clearvars -except constt box Rmm real_gr compl_gr period RSQ 
    % EDIT THE FILE NAME DIRECTORY
    box{constt} = 32;
    filename_in  = ['/data/novadisk/vs391/snoopy_kinematic_dynamo/kd_box_' num2str(box{constt}) '/snoopy/kinematicOutput/timevar' num2str(constt) '.dat'];
    filename_out = ['/data/novadisk/vs391/snoopy_kinematic_dynamo/kd_box_' num2str(box{constt}) '/snoopy/kinematicOutput/tim_sum' num2str(box{constt}) '.dat'];
    
    try 
        % IMPORT THE DATA TO READ IN COLUMN BY COLUMN 
        full_file=importdata(filename_in);
        timevar=transpose(full_file.data);
        nvar=size(timevar,1);
        var_name=strread(full_file.textdata{1},'%s',nvar);

        % ASSIGN VARIABLE WITH HEADINGS 
        for ii=1:nvar
            assignin('base',var_name{ii},timevar(ii,:)); 
        end

        if exist('em') && exist('Rm') && exist('t')
            
            % @em_tr    GET RID OF THE INITIAL DIFFUSION 
            % @ln_em_tr TAKE LN OF EM
            if (size(em(em>1.0))>1)
                em_tr = em(em>1.0);
                ln_em = log(em_tr);
                tt    = t(em>1.0);
            else 
                em_tr = em(t>0);
                ln_em = log(em_tr);
                tt    = t(t>0);
            end
            % LINEAR FIT TO THE LN(EM)
            % @f(1) IS THE REAL EXP GROWTH RATE 
            f = polyfit(tt,ln_em,1);
            real_gr{constt} = f(1);  
        
            % @RSQ - AVERAGE SUM OF RESIDUALS SQUARED
            % IF @RSQ -> 0 -- GOOD LINEAR FIT
            yfit = polyval(f,tt);
            yresid = ln_em - yfit;
            if (size(tt,2) > 2) 
                RSQ{constt} = sum(yresid.^2)/(size(tt,2)-2);   
            else 
                RSQ{constt} = 0;
            end
            
            if (RSQ{constt} > 0.2)
                % GET THE OSCILLITARY PART BY EXP(-f(1)*tt)
                em_osc = em_tr .* exp(-f(1) .* tt);

                % AUTOCORRELATION
                % NORMALISE AND SUBTRACT THE MEAN
                em_osc_m = em_osc - mean(em_osc);
                [autocor,lags] = xcorr(em_osc_m);
                % SINCE AUTOcorrelation, LOOK AT THE POSITIVE LAG
                autocor = autocor(lags>0)./max(autocor);
                lags = lags(lags>0);
                % FIND EXTREMA POINTS
                if (size(autocor)<2)
                    period{constt}   = 0.0;
                    compl_gr{constt} = 0.0;
                else
                    dr = diff(autocor);
                    eIdx = find(dr(1:end-1) .* dr(2:end) <= 0) + 1;

                    % CHECK THE FUNCTION HAS EXTREMA
                    if (isempty(eIdx)==0)
                        [~,loc] = sort(autocor(eIdx), 'descend');
                        loc = loc(1:min(2,end));  
                        % convert lag to real time 
                        lag_step = (tt(end)-tt(1))/(size(tt,2)-1);
                        % the period is given by the max converted lag
                        period{constt} = lags( eIdx(loc(1)))*lag_step;
                        compl_gr{constt}=2*pi/period;
                    else 
                        period{constt}   = 0.0;
                        compl_gr{constt} = 0.0;
                    end
                end
            else
                period{constt}   = 0.0;
                compl_gr{constt} = 0.0;
            end
           
        % SAVE RM ENTRIES
        fprintf('In timevar%d: Rm=%f\n',constt, Rm(1));
        Rmm{constt} = Rm(1);
        end
    catch ME
        fprintf('Error with timevar%d.dat file %s\n',constt, ME.message);
        continue; 
    end 

   
end
% Save the results to the table   

tblC_title={'Box'; 'Rm'; 'RGR'; 'RSq_r';'CGR';'Period'};
tblC = table(transpose(box),transpose(Rmm),transpose(real_gr),transpose(RSQ),transpose(compl_gr),transpose(period),'VariableNames', transpose(tblC_title));
writetable(tblC, filename_out,'Delimiter', 'tab');



