galloway1=[ 8 -0.0076    0.628   10];
galloway2=[ 10 0.0056     0.628   10 ];
galloway3=[ 16 0.0009     0.6     10.5];
galloway4=[ 20 -0.0032    0.565   11.1];
galloway5=[ 24 -0.023     0.548   11.44];
galloway6=[ 27 0          0.128   49];
galloway7=[ 30 0.013      0.124   50.6];
galloway8=[ 50 0.040      0.097   64.5];

% Edit next two lines 
galloway=[galloway3];
full_file=importdata(['timevar3.dat']);

clear galloway1; clear galloway2;clear galloway3;clear galloway4;
clear galloway5; clear galloway6;clear galloway7;clear galloway8;

timevar=transpose(full_file.data);
nvar=size(timevar,1);
var_name=strread(full_file.textdata{1},'%s',nvar);

for i=1:nvar
   assignin('base',var_name{i},timevar(i,:)); 
end


% @em_tr    GET RID OF THE INITIAL DISPERSION 
% @ln_em_tr TAKE LN OF EM
em_tr = em(17:end);
ln_em = log(em_tr);
tt = t(17:end);

% LINEAR FIT TO THE LN(EM)
% @f(1) IS THE REAL EXP GROWTH RATE 
f = polyfit(tt,ln_em,1);
yfit = polyval(f,tt);
yresid = ln_em - yfit;
sum(yresid.^2)/size(tt,2)

% GET THE OSCILLITARY PART BY EXP(-f(1)TT)
em_osc = em_tr .* exp(-f(1) .* tt);


% AUTOCORRELATION
% Normalise by subtracting the mean 
em_osc_m = em_osc - mean(em_osc);
[autocor,lags] = xcorr(em_osc_m);
% since AUTOcorrelation, look at positive lag
autocor = autocor(lags>0)./max(autocor);
lags = lags(lags>0);
% find extrema points
dr = diff(autocor);
eIdx = find(dr(1:end-1) .* dr(2:end) <= 0) + 1;
[~,loc] = sort(autocor(eIdx), 'descend');
loc = loc(1:min(3,end));  
% convert lag to real time 
lag_step = (tt(end)-tt(1))/(size(tt,2)-1);
% the period is given by the max converted lag
period = lags( eIdx(loc(3)))*lag_step;
compl_gr=2*pi/period;
       		
% PLOT

figure('units','normalized','position',[.1 .1 .5 .4])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',14) 

subplot(2,3,1);
plot(t, log(em));
xlabel('time'); title('ln(EM)'); 
xlim([0,200]);

subplot(2,3,2); 
plot(tt, ln_em, 'b',tt,f(1).*tt+f(2), 'r');
xlabel('time');title('trancated ln(EM)');
legend('ln(EM)', 'Linear fit');
xlim([0,200]);

subplot(2,3,3);
plot(tt, em_osc);
xlabel('time'); title('osc EM');
xlim([0,200]);

subplot(2,3,[4,5]);
plot(lags.*lag_step,autocor);
xlabel('Time'); title('Autocorrelation');


% TABLE
format long
sdata=[Rm(1) galloway(1); f(1)/2. galloway(2); compl_gr/2. galloway(3); period*2. galloway(4)];

rnames={'Rm','Real GR', 'Complex GR', 'Period'};
colnames= {'SNOOPY','GALLOWAY'};
t = uitable('Data', sdata, 'ColumnName', colnames,'RowName',rnames,'Position',[650 97 1 1] );
tableextent = get(t,'Extent');
oldposition = get(t,'Position');
newposition = [oldposition(1) oldposition(2) tableextent(3) tableextent(4)];
set(t, 'Position', newposition);
