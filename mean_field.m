for const=1:10
try 
    filename_out = ['/Users/valeriashumaylova/Dropbox/matlab/box_5/data/positive_' num2str(const)];
    filename_in  = ['/Users/valeriashumaylova/Dropbox/matlab/box_5/data/modes_' num2str(const)];

    full_file=importdata(filename_in);

    timevar=full_file.data;
    nvar=size(timevar,2);

    % Create a table with all variables
    tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6),timevar(:,7));
    KX=timevar(:,1);
    KY=timevar(:,2);
    KZ=timevar(:,3);
    RG=timevar(:,4);
    CG=timevar(:,6);

    tblB_title={'KX'; 'KY'; 'KZ';'RGR'; 'RSq_r';'CGR';'Period'};
    for k=1:size(KX)
        if (timevar(k,4) > 0.0) 
            tblB{k,:} = table(KX(k), KY(k), KZ(k), timevar(k,4), timevar(k,5),timevar(k,6),timevar(k,7), 'VariableNames', transpose(tblB_title));
        end 
    end 

    writetable(tblB, filename_out,'Delimiter', 'tab');
catch ME
        fprintf('Error modes_%d file %s\n', const, ME.message);
        continue; 
    end 
end

        