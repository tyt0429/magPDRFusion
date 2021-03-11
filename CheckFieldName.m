function [ds] = CheckFieldName(ds)
for i=1:length(ds.VariableNames)
    if strfind(ds.VariableNames{i},'GRAVITY')
        if strfind(ds.VariableNames{i},'GRAVITYX')
            ds.VariableNames{i}='Gra_x';
        elseif strfind(ds.VariableNames{i},'GRAVITYY')
            ds.VariableNames{i}='Gra_y';
        elseif strfind(ds.VariableNames{i},'GRAVITYZ')
            ds.VariableNames{i}='Gra_z';
        end
    elseif strfind(ds.VariableNames{i},'MAGNETIC')
        if strfind(ds.VariableNames{i},'X')
            ds.VariableNames{i}='Mag_x';
        elseif strfind(ds.VariableNames{i},'Y')
            ds.VariableNames{i}='Mag_y';
        elseif strfind(ds.VariableNames{i},'Z')
            ds.VariableNames{i}='Mag_z';
        end
        
    elseif strfind(ds.VariableNames{i},'ACCELEROMETER')
        if strfind(ds.VariableNames{i},'X')
            ds.VariableNames{i}='Acc_x';
        elseif strfind(ds.VariableNames{i},'Y')
            ds.VariableNames{i}='Acc_y';
        elseif strfind(ds.VariableNames{i},'Z')
            ds.VariableNames{i}='Acc_z';
        end
        
    elseif strfind(ds.VariableNames{i},'GYROSCOPE')
        if strfind(ds.VariableNames{i},'GYROSCOPEX')
            ds.VariableNames{i}='Gyr_x';
        elseif strfind(ds.VariableNames{i},'GYROSCOPEY')
            ds.VariableNames{i}='Gyr_y';
        elseif strfind(ds.VariableNames{i},'GYROSCOPEZ')
            ds.VariableNames{i}='Gyr_z';
        end
    elseif strfind(ds.VariableNames{i},'ORIENTATIONZ')
        ds.VariableNames{i}='Ori_z';
        
    end
end

