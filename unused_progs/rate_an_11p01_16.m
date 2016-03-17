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

clear all;
%full_file=importdata('rates_16_6/rates_r00.txt');
full_file=importdata('../../rates_r00.txt');

timevar=full_file.data;
nvar=size(timevar,2);
%var_name=strread(full_file.textdata{2},'%s',nvar);

% Assign variables to columns of data
%for i=1:nvar
   %assignin('base',var_name{i},timevar(:,i));
%end

% Create a table with all variables
%tblA = table(KX,KY,KZ,TIME,EM,GR);

tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
KX=timevar(:,1);
% Sort the rows of the table
tblB = sortrows(tblA);

% Find how many different modes are present
% @mode_box - array that contains the location of last mode of the group
% @n_modes  - number of various dynamo modes
i=0;
k=0;
mode_box{1}=1;
for j=1:(size(KX,1)-1)
    if (tblB{j,1}==tblB{j+1,1} &tblB{j,2}==tblB{j+1,2} & tblB{j,3}==tblB{j+1,3})
        i=i+1;
    else
        k = k+1;
        mode_box{k+1}=j+1;
    end
end
mode_box{k+2}=size(KX,1);
n_modes = size(KX,1) - i;

% Look at each mode and determine whether its growing exponentially
% and if it is exserting an oscillating bahaviour
for i=1:n_modes
    clear x;            clear y;        clear x_arr;        clear y_arr;
    clear y_arr_osc_sh; clear yfit;     clear yresid;       clear residuals;
    clear ln_y;         clear ln_y_arr; clear x_0;          clear x_arr;
    clear x0;           clear x00;      clear y_0;          clear y_arr;
    clear x_arr_sh;     clear y_arr_sh; clear ln_y_arr_sh;  clear objfcn;
    % Check that a mode has data at more than one time
    % for such modes, take the growth rate from Snoopy
   if (mode_box{i+1}-mode_box{i} == 1)
            gr_arr{i}=0.0;
            RSq{i}=0.0;
            lin_gr{i}=0.0;
            lin_Rsq{i}=0.0;
            cgr_arr{i}=0.0;
            cplx_RSq{i}=0.0;
   else
       % Get the EM and time data to fit an exp
       % @x is time
       % @y is B(k)
       for j=1:(mode_box{i+1}-mode_box{i})
           x{j} = tblB{j+mode_box{i}-1,4};
           y{j} = tblB{j+mode_box{i}-1,5};
       end


       % @ln_y is ln(B(k))
       % @x_0  is time (ommitting the zero B(k))
       % @y_0  is B(k) (ommitting the zero B(k))
       % @kk   is a dummy variable
       kk=0;
       for k=1:size(y,2)
            if (y{k} > 0.0 & isfinite(y{k})==1)
                kk=kk+1;
                ln_y{kk} = log(y{k});
                x_0{kk} = x{k};
                y_0{kk} = y{k};
            end
       end
       clear kk;

       % Convert tables x,y and ln_y into arrays
       ln_y_arr = cell2mat(ln_y);
       x_arr = cell2mat(x_0);
       y_arr = cell2mat(y_0);

       % not sure x_arr_sh(ii)    = x_arr(kk);

       % Allow the simulation to run for some time and look at the slope
       % for times (tf/2, tf), where tf is the final time.
       % @x_arr_sh time for the 2nd half of simulation
       % @y_arr_sh B(k) for the 2nd half of simulation
       % @ln_y_arr_sh ln(B(k)) for the 2nd half of simulation
       if (size(x_arr,2)>5)
           ii=0;
           if (mod(size(x_arr,2),2)== 0)
               for kk= int64((size(x_arr,2)/2))  :    size(x_arr,2)
                    ii=ii+1;
                    x_arr_sh(ii)    = x_arr(kk);
                    y_arr_sh(ii)    = y_arr(kk);
                    ln_y_arr_sh(ii) = ln_y_arr(kk);
               end
           else
               for kk= int64((size(x_arr,2)-1)/2):size(x_arr,2)
                    ii=ii+1;
                    x_arr_sh(ii)    = x_arr(kk);
                    y_arr_sh(ii)    = y_arr(kk);
                    ln_y_arr_sh(ii) = ln_y_arr(kk);
               end
           end
       else
           x_arr_sh     = x_arr;
           y_arr_sh     = y_arr;
           ln_y_arr_sh  = ln_y_arr;
       end
       clear ii;

       if (size(x_arr_sh) < 3)
            gr_arr{i}=0.0;
            RSq{i}=0.0;
            lin_gr{i}=0.0;
            lin_Rsq{i}=0.0;
            correl{i}=0.0;
            period{i}=0.0;
       else
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
       RSq{i} = sum(yresid.^2);

       clear yresid;

	% If the linear fit is good, then
	% ignore next steps. Otherwise,
	% calculate the complex growth rate
    
    % Multiply EM by exp(-real gr * time) to get a non-growing
    % oscillating solution to find the complex gr rate.
    % @y_arr_osc_sh is the oscillating solution B(k)exp(-real gr t)
    for jj=1:size(y_arr_sh,2) 
        y_arr_osc_sh(jj) = y_arr_sh(1,jj).*exp(-f(1)*x_arr_sh(1,jj));
    end
            
        if (RSq{i} < 0.01 | max(y_arr_osc_sh) < 1e-5*gr_arr{i}) 
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
       		y_arr_osc_sh = y_arr_osc_sh - mean(y_arr_osc_sh);
            % Calculate the smallest time interval in your data
            for jj=1:(size(x_arr_sh,2)-1)
                t_int(jj) = x_arr_sh(jj+1) - x_arr_sh(jj);
            end
            t_min = min(t_int);
            clear t_int; clear jj;
            % Fill in the missing gaps with normalised @y_arr_osc_sh = 0
            % This will not make a difference in correlation
            y_arr_osc_sh_0(1:round(x_arr_sh(end))/round(t_min)) = 0.0; 
            for jj = 1:round(x_arr_sh(end))/round(t_min)        %1:1298/2=649
                time(jj)=jj*round(t_min);
                for kk =1:size(x_arr_sh,2)                      %1:104        
                    if (jj*round(t_min) ==  round(x_arr_sh(kk)))             %870:1298
                        %disp(['jj=',num2str(jj),' x=',num2str(x_arr_sh(kk)), ' y=',num2str(y_arr_osc_sh(kk))])
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
            autocor = autocor(lags>=0)./max(autocor);
            lags = lags(lags>=0);

		    % find extrema points
		    dr = diff(autocor);
		    eIdx = find(dr(1:end-1) .* dr(2:end) <= 0) + 1;
		    [~,loc] = sort(autocor(eIdx), 'descend');
		    %loc = loc(1:min(2,size(loc)));                     %# take the highest 2 values

            % Convert lags to real time
            lag_step = round(t_min);
            % the period is given by the max converted lag
            period{i} = lags( eIdx(loc(2)))*lag_step;
            % Complex growth rate is then given by
       		correl{i}=2*pi/period{i};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       		%%%%%%%%%%%% PLOT %%%%%%%%%%%%%
       		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       		figure(i)
       		subplot(3,1,1)
       		plot(time(round(x_arr_sh(1)/t_min):end), y_arr_osc_sh_0(round(x_arr_sh(1)/t_min):end),'r--*');
     	  	title('Oscillating B(t)');

            subplot(3,1,2)
       		autocor = autocor(lags>=0)./max(autocor);
       		lags = lags(lags>=0).*round(t_min);
       		plot(lags,autocor);
       		xlabel('Time');
       		ylabel('Autocorrelation');
            
       		%subplot(3,1,3)
       		%plot(time(round(x_arr_sh(1)/t_min):end),y_arr_osc_sh_0(round(x_arr_sh(1)/t_min):end),'r',time(round(x_arr_sh(1)/t_min):end), max(y_arr_osc_sh_0).*sin(2*pi/correl{i}*time(round(x_arr_sh(1)/t_min):end)) );
      	 	%ylabel('sin of autocorr');
       		
    end
    end
	
       		
            
	end
   end
  
% Save the results to the table
for i=1:n_modes
   kx_arr{i} = tblB{mode_box{i},1};
   ky_arr{i} = tblB{mode_box{i},2};
   kz_arr{i} = tblB{mode_box{i},3};
end

tblC_title={'KX'; 'KY'; 'KZ';'RGR'; 'RSq_r';'CGR';'Period'};
tblC = table(transpose(kx_arr),transpose(ky_arr),transpose(kz_arr),transpose(gr_arr),transpose(RSq),transpose(correl),transpose(period),'VariableNames', transpose(tblC_title));
tblD = sortrows(tblC,-4);
