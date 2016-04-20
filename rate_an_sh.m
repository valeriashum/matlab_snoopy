%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	This file is part of the Snoopy code.

%   Snoopy code is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    Snoopy code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the data
% @i is only used for outputting to the table
% save this file in KinematicOutput/graphics
rep='/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/snoopy1/kinematicOutput';
rep1='/media/vs391/TOSHIBAHARD/WORK/u_abc/box_10_kinematicOutput_5/';
 for constt = 2  % file number
    for const = 0:7    %  number of cores
        fprintf('Begin: rates%d/rates_r0%d\n', constt, const);
         clearvars -except const constt folder rep rep1
         box_ar = 10;
         filename_out = [rep1, '/rates' num2str(constt) '/modes_r0' num2str(const) '.txt'];
         filename_in  = [rep1, '/rates' num2str(constt) '/rates_r0' num2str(const) '.txt'];
         filename_rm  = [rep1, '/timevar' num2str(constt) '.dat'];

        try
            % get the magnetic Reynolds number for the sim
            rm_file=importdata(filename_rm);
            rmvar=rm_file.data;
            nrm=size(rmvar,2);

            % read in from rates0*        
            full_file=importdata(filename_in);

            timevar=full_file.data;
            nvar=size(timevar,2);

            % Create a table with all variables
            %tblA = table(KX,KY,KZ,TIME,EM,GR);

            tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
            KX=timevar(:,1);
            % Sort the rows of the table
            tblB = sortrows(tblA);

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
            
            for j=1:(size(KX,1)-1)
                if (tblB{j,1}==tblB{j+1,1} & tblB{j,2}==tblB{j+1,2} & tblB{j,3}==tblB{j+1,3})
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
            mode_box{k+2}=size(KX,1);

            % Calculate total number of modes growing in the system
            % Should be equal to size(mode_box,2)-1
            n_modes = size(KX,1) - iii -n_sin;
            
            % Preallocate the arrays you submit to the table in the end
            gr_arr = cell(1,n_modes);
            RSq= cell(1,n_modes);
            correl= cell(1,n_modes);
            period= cell(1,n_modes);
            % Look at each mode and determine whether its growing exponentially
            % and if it is exserting an oscillating bahaviour
            fprintf('UPDATE: Document has been uploaded\n');   
            clear j k n_sin iii;
            
            % Maximum power at a fourier mode throughtout the simulation is
            % given by  
            
            for j=1:(size(KX,1))
                em_cell{j} = tblB{j,5};
                time_cell{j} = tblB{j,4};
            end
            em = cell2mat(em_cell);
            time = cell2mat(time_cell);
            max_power = max(em); 
            max_time = max(time);
            % Set power limit:
            % (power-max)/max < limit is good  
            power_limit = 0.5;
            
            % Set a constraint on life-span of modes investigated: 
            % If max_time - time < time_limit -> good enough
            time_limit = 0.;
            % Set a constraint on mean sum of residuals squared: 
            % If R^2 < res_limit -> good enough
            res_limit=1e-6;
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP      %%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:n_modes
                clear x;            clear y;        clear x_arr;        clear y_arr;
                clear y_arr_osc_sh; clear yfit;     clear yresid;       clear residuals;
                clear ln_y;         clear ln_y_arr; clear x_0;          clear x_arr;
                clear x0;           clear x00;      clear y_0;          clear y_arr;
                clear x_arr_sh;     clear y_arr_sh; clear ln_y_arr_sh;  clear objfcn;
                % Check that a mode has data at more than one time
                % for such modes, take the growth rate from Snoopy

                % Set a min number of points to fit a line through
                min_number = 4;

               if (mode_box{i+1}-mode_box{i} == 1)
                        gr_arr{i}=0.0;
                        RSq{i}=0.0;
                        correl{i}=0.0;
                        period{i}=0.0;
               else
                   % Get the EM and time data to fit an exp
                   % @x is time
                   % @y is B(k)
                   for j=1:(mode_box{i+1}-mode_box{i})
                       x{j} = tblB{j+mode_box{i}-1,4};
                       y{j} = tblB{j+mode_box{i}-1,5};
                   end
                   clear j;
                   
                   % Get rid of repeating entries (rounding)
                   % Since @x(1) = 0, which we eliminate later,
                   % we can start with k>1 and set all unwanted 
                   % entries to zero.
                   for j= 2:(mode_box{i+1}-mode_box{i})
                       if (round(x{j}) == round(x{j-1}))
                           y{j} = 0.0;
                       end
                   end
                   clear j;
                   
                   % @ln_y is ln(B(k))
                   % @x_0  is time (ommitting the zero B(k))
                   % @y_0  is B(k) (ommitting the zero B(k))
                   % @kk   is a dummy variable
                   
                   jj=0;
                   for j= 1:(mode_box{i+1}-mode_box{i})
                        if (y{j} > 0.0 && isfinite(y{j})==1)
                            jj=jj+1;
                            ln_y{jj} = log(y{j});  
                            x_0{jj} = x{j};
                            y_0{jj} = y{j};
                        end
                   end
                   clear j jj;


                   % Convert tables x,y and ln_y into arrays
                   ln_y_arr = cell2mat(ln_y);
                   x_arr = cell2mat(x_0);
                   y_arr = cell2mat(y_0);
                   
                   % Only find growth rates for modes that actually grow
                   % in the long run 
                   
                   %if ((y_arr(end)-max_power)/max_power < power_limit | max_time-x_arr(end) > time_limit )
                   if (size(x_arr,2)< 2)
                       gr_arr{i}=0.0;
                        RSq{i}=0.0;
                        Bntfin{i} = y_arr(end);
                        period{i} = 0.0;
                        correl{i} = 0.0;
                   else
                   if (y_arr(end) < y_arr(end-1)| max_time-x_arr(end) > time_limit )
                       gr_arr{i}=0.0;
                       RSq{i}=0.0;
                       Bntfin{i} = y_arr(end);
                       period{i} = 0.0;
                       correl{i} = 0.0;
                       %fprintf('For kx=%f, ky=%f, kz=%f,y_arr(2)=%e,y_arr(end)=%e\n', tblB{mode_box{i},1},tblB{mode_box{i},2},tblB{mode_box{i},3},y_arr(2), y_arr(end));
                   else
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%  REAL GR   %%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       % If there are less than 3 points, we cannot 
                       % calculate any of the growth rates 
                       if (size(x_arr,2)< min_number)
                            gr_arr{i}=0.0;
                            RSq{i}=0.0;
                            Bntfin{i} = y_arr(end);
                            period{i} = 0.0;
                            correl{i} = 0.0;
                       elseif (size(x_arr,2)== min_number)

                            x_arr_sh = x_arr;
                            y_arr_sh = y_arr;
                            ln_y_arr_sh = ln_y_arr;     

                            % Linear fit to TIME vs ln(EM), real growth rate is given
                            % by the slope of the line of best fit.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%% POLYFIT %%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % p(x)=p1xn+p2xn−1+...+pnx+pn+1.
                            % @f(1) the slope on line fit for n=1
                            f = polyfit(x_arr_sh,ln_y_arr_sh,1);
                            gr_arr{i} = double(f(1));

                            % Find how good of a fit the linear curve is
                            % by looking at R^2, which predicts % variance in EM
                            % For a good linear fit:
                            % @RSq     -> 0
                            yfit = polyval(f,x_arr_sh);
                            yresid = ln_y_arr_sh - yfit;
                            RSq{i} = sum(yresid.^2)/(size(x_arr_sh,2)-2);
                            Bntfin{i} = y_arr_sh(end);
                       else 

                            % Since the time evolution for each fourier component of
                            % magnetic field is different, look at various sets of data
                            % points to find the one that gives the lowest sum of residuals
                            % squared: 
                            % @jjj is the number of points - 3 for each Bn, since 3 is the minimum
                            % number of data points that we can fit a line through
                            % @x_arr_sh - time 
                            % @y_arr_sh -  B(k) 
                            % @ln_y_arr_sh - ln(B(k))
                            for jjj = 1:size(x_arr,2)-min_number 

                                clear x_arr_sh; clear y_arr_sh;  clear ln_y_arr_sh; 

                                ii=0;
                                for kk = jjj:size(x_arr,2)
                                    ii=ii+1;
                                    x_arr_sh(ii)    = x_arr(kk);
                                    y_arr_sh(ii)    = y_arr(kk);
                                    ln_y_arr_sh(ii) = ln_y_arr(kk);
                                end
                                clear ii kk;    

                                % Linear fit to TIME vs ln(EM), real growth rate is given
                                % by the slope of the line of best fit.
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%% POLYFIT %%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % p(x)=p1xn+p2xn−1+...+pnx+pn+1.
                                % @f(1) the slope on line fit for n=1
                                f = polyfit(x_arr_sh,ln_y_arr_sh,1); 
                                lin_gr(jjj) = double(f(1));

                                % Find how good of a fit the linear curve is
                                % by looking at R^2, which predicts % variance in EM
                                % For a good linear fit:
                                % @RSq     -> 0
                                yfit = polyval(f,x_arr_sh);
                                yresid = ln_y_arr_sh - yfit;
                                lin_resid(jjj) = sum(yresid.^2)/(size(x_arr_sh,2)-2);

                                if (lin_resid(jjj) < res_limit)
                                    break
                                end        
                            end

                            % Find which time data points give us the best linear fit
                            % by determining the min of R^2_residuals
                            % @ind_min_res is the index of the min sum of residuals but
                            
                            [min_res,ind_min_res] = min(lin_resid);
                            RSq{i} = min_res;
                            gr_arr{i} = lin_gr(ind_min_res);
                            
                            % @Bntfin - power in a fourier component at
                            % final time from the rates file
                            Bntfin{i} = y_arr_sh(end);
                            
                            % @jjj is used to redefine @x_arr_sh etc
                            % that were used to find the linear growth rate
                            ii=0;
                            for kk = ind_min_res:size(x_arr,2)
                                ii=ii+1;
                                x_arr_sh(ii)    = x_arr(kk);
                                y_arr_sh(ii)    = y_arr(kk);
                                ln_y_arr_sh(ii) = ln_y_arr(kk);
                            end
                            clear ii kk;

                            %fprintf('For kx=%f\tky=%f\tkz=%f:\t%d:%d, \t gr=%e\n',tblB{mode_box{i},1},tblB{mode_box{i},2},tblB{mode_box{i},3},ind_min_res,size(x_arr,2),gr_arr{i} );
                            %fprintf('For kx=%f\tky=%f\tkz=%f:\t%d:%d, \t gr=%e\n',tblB{mode_box{i},1},tblB{mode_box{i},2},tblB{mode_box{i},3},round(x_arr(ind_min_res)),round(x_arr(end)),gr_arr{i} );
                            clear lin_resid;  clear lin_gr;  clear min_res;   clear  ind_min_res; 
                       end
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLEX GR %%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       if (size(x_arr,2)< min_number)
                            correl{i}=0.0;
                            period{i}=0.0;
                       else       
                            % Multiply EM by exp(-real gr * time) to get a non-growing
                            % oscillating solution to find the complex gr rate.
                            % @y_arr_osc_sh is the oscillating solution B(k)exp(-real gr t)
                            for jj=1:size(y_arr_sh,2) 
                                y_arr_osc_sh(jj) = y_arr_sh(jj).*exp(-gr_arr{i}*x_arr_sh(jj));
                            end
                            clear jj;

                            % If the linear fit is good, then
                            % ignore next steps. Otherwise,
                            % calculate the complex growth rate
                            if (RSq{i} < 0.01 | max(y_arr_osc_sh) < 1e-10*gr_arr{i}) 
                                correl{i}=0.0;
                                period{i}=0.0;
                            else
                                % Try to find the period/complex growth rate
                                % using the cross-correlation
                                % http://uk.mathworks.com/help/signal/ug/find-periodicity-using-autocorrelation.html
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%  AUTOCORRELATION %%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % Normalise by subtracting the mean
                                y_arr_osc_sh(1:end) = y_arr_osc_sh(1:end) - mean(y_arr_osc_sh);
                                % Calculate the smallest time interval in your data
                                for jj=1:(size(x_arr_sh,2)-1)
                                    t_int(jj) = x_arr_sh(jj+1) - x_arr_sh(jj);
                                end
                                clear jj; 

                                % Because of a glitch that appears due to the output time
                                % change, SNOOPY outputs data at e.g. 4.0, 4.0001 and then 
                                % switches to a new t_out. To avoid calculating the autocor
                                % for each e.g. 0.0001, we ignore the values at 4.0001. 
                                % In order to do so, we need to alter the min t_out. 
                                if (round(min(t_int)) == 0)
                                    [t_int_s,loc] = sort(t_int(t_int>0.51), 'ascend');
                                    t_min = t_int_s(1);
                                    clear t_int_s;
                                else 
                                    t_min = min(t_int);
                                end 
                                clear t_int; clear jj;

                                % Fill in the missing gaps with normalised @y_arr_osc_sh = 0
                                % This will not make a difference in correlation
                                y_arr_osc_sh_0(1:round(x_arr_sh(end)/t_min)) = 0.0; 
                                for jj = 1:round(x_arr_sh(end)/t_min)       
                                    time(jj)=jj*round(t_min);
                                    for kk =1:size(x_arr_sh,2)                             
                                        if (jj ==  round(x_arr_sh(kk)/t_min))             
                                            y_arr_osc_sh_0(jj) = y_arr_osc_sh(kk);   
                                        end
                                    end
                                end

                                % Calculate the autocorrelation for  @y_arr_osc_sh_0
                                % Lags can be converted into time by
                                % time = @lags * @lag_step
                                [autocor,lags] = xcorr(y_arr_osc_sh_0);

                                % Since we are finding the AUTOcorrelation, look at the positive lags,
                                % which is a mirror reflection of negative lags
                                % Normalise the correlation
                                autocor = autocor(lags>0)./max(autocor);
                                lags = lags(lags>0);

                                % find extrema points
                                dr = diff(autocor);
                                eIdx = find(dr(1:end-1) .* dr(2:end) <= 0) + 1;

                                if (isempty(eIdx)==0)
                                    [~,loc] = sort(autocor(eIdx), 'descend');
                                    % Convert lags to real time
                                    lag_step = round(t_min);
                                    % The period is given by the max converted lag 
                                    % Take the highest 1 value since the zero entry of autocor is 0
                                    period{i} = lags( eIdx(loc(min(1,size(loc,2)))))*lag_step;     
                                    % Complex growth rate is then given by
                                    correl{i}=2*pi/period{i};

                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%%%%%%%%%%% PLOT %%%%%%%%%%%%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    %figure(i)
                                    %subplot(3,1,1)
                                    %plot(time(round(x_arr_sh(1)/t_min):end), y_arr_osc_sh_0(round(x_arr_sh(1)/t_min):end),'r--*');
                                    %title('Oscillating B(t)');

                                    %subplot(3,1,2)
                                    %autocor = autocor(lags>=0)./max(autocor);
                                    %lags = lags(lags>=0).*round(t_min);
                                    %plot(lags,autocor);
                                    %xlabel('Time');
                                    %ylabel('Autocorrelation'); 

                                    %subplot(3,1,3)
                                    %plot(time(round(x_arr_sh(1)/t_min):end),y_arr_osc_sh_0(round(x_arr_sh(1)/t_min):end),'r',time(round(x_arr_sh(1)/t_min):end), max(y_arr_osc_sh_0).*sin(2*pi/correl{i}*time(round(x_arr_sh(1)/t_min):end)) );
                                    %ylabel('sin of autocorr');
                                else
                                    period{i} = 0.0;
                                    correl{i} = 0.0;
                                end
                            end
                       end
                   end 
                   end
               end
                % Save the results to the table   

                kx_arr{i} = tblB{mode_box{i},1};
                ky_arr{i} = tblB{mode_box{i},2};
                kz_arr{i} = tblB{mode_box{i},3};
                rm{i} = rmvar(1,3);
                box{i}= box_ar;
            end
            
        % Normalise fourier coeffs B_n by max(B_n)    
        bn = cell2mat(Bntfin); 
        max_bn= max(bn);
        clear bn;
        
        % Save to the table 
        tblC_title={'Box';'Rm';'KX'; 'KY'; 'KZ';'RGR'; 'RSq_r';'Bntfinnorm';'CGR';'Period'};
        tblC = table(transpose(box),transpose(rm), transpose(kx_arr),transpose(ky_arr),transpose(kz_arr),transpose(gr_arr),transpose(RSq),transpose(Bntfin),transpose(correl),transpose(period),'VariableNames', transpose(tblC_title));
        tblD = sortrows(tblC,-6);
        
%         % To analyse the data in @tblD, first delete the rows with zero 
%         % growth rate and normalise fourier coeffs B_n by max(B_n)
%         for i=1:n_modes
%             if (tblD.RGR{i} <= 0.0) 
%                 tblD(i,:)=[];
%             end 
%         end
        
        
        % OUTPUT
        writetable(tblD, filename_out,'Delimiter', 'tab');

        catch ME
            fprintf('Error with the %d rates %d file %s\n',constt, const, ME.message);
            continue; 
        end 
        
        
    
    end
 end
