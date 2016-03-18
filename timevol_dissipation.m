
files1=[    {'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/results/nonlin_box1_rm12_re12_implicit/timevar'}... 
            {'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/results/nonlin_box1_rm12_re12_implicit_res32/timevar'}... 
            {'/store/ASTRO/vs391/nonlinear_dynamo/u_abc/results/nonlin_box1_rm12_re12/timevar'}      ];
legendInfo1 = [     {'$\textit{Implicit Dissipation } 16^3$'}...
                    {'$\textit{Implicit Dissipation } 32^3$'}...
                    {'$\textit{Explicit Dissipation } 32^3$'}];    
                
sizearray=[size(files1,2) ];    
maxsize=max(sizearray);
         
%for co1=1:maxsize
 %   color{co1}=rand(1,3); 
%end
color{1}=[0 0 1]
color{2}=[1 0 0]
color{3}=[0 0 0]
for j=1:size(files1,2)
    full_file=importdata(files1{1,j});
    timevar=transpose(full_file.data);
    nvar=size(timevar,1);
    var_name=strread(full_file.textdata{1},'%s',nvar);

    for i=1:nvar
       assignin('base',var_name{i},timevar(i,:)); 
    end
    figure(1)
    if exist('t')&&exist('ev')&&exist('em')
        subplot(2,1,1)
            plot(t,ev,'color',color{j},'LineWidth',1.5);
            hold on;
            title('$\textit{ABC flow at } R_m = R_e=12$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_K$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo1,'Location','northwest','Interpreter', 'latex','fontsize',20)
            xlim([0 3000])
        subplot(2,1,2)
            plot(t,em,'color',color{j},'LineWidth',1.5);
            hold on;
            title('$\textit{ABC flow at } R_m = R_e=12$','fontsize',16, 'Interpreter', 'latex');
            xlabel('$T$','fontsize',16, 'Interpreter', 'latex');
            ylabel('$E_M$','fontsize',16, 'Interpreter', 'latex');
            legend(legendInfo1,'Location','northwest','Interpreter', 'latex','fontsize',20)
            xlim([0 3000])
    end            
end 

