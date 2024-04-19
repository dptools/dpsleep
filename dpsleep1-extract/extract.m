function extract(read_dir)

% Check if the path is properly formatted
if ~ endsWith(read_dir, '/')
    read_dir = strcat(read_dir, '/');
end

% Get the files from directory
files = dir(strcat(read_dir, '*.csv'));
files_len = length(files);

% Exit if there are no files to read
if files_len == 0
    display('Files do not exist under this directory.');
    exit(1);
end

% For every csv file under the directory, create a .mat file
for i=1:files_len
    file_path = strcat(files(i).folder, '/', files(i).name);
    file_name = erase(files(i).name, '.csv');

    % Try reading in the file
    try
        % Skipping the csv header and timestamp column
        dt = csvread(file_path, 1, 1);
        dt2 = dt(:, 1:end-1);
        fs = dt(:, end);
        dtt1=dt2(:,1);
        dtt2=dt2(:,2);
        dtt3=dt2(:,3);
        dtt4=dt2(:,4);
        dtt5=dt2(:,5);
        dtt6=dt2(:,6);
        ldt=length(dtt1);
        dtv1=reshape(dtt1,ldt/24,24);
        dtv2=reshape(dtt2,ldt/24,24);
        dtv3=reshape(dtt3,ldt/24,24);
        dtv4=reshape(dtt4,ldt/24,24);
        dtv5=reshape(dtt5,ldt/24,24);
        dtv6=reshape(dtt6,ldt/24,24);
        for dy=2:24
            dty=dtv1(:,dy);
            dty0=dtv1(:,dy-1);
            if (~isnan(dty(1)) && isnan(dty(end))) && isnan(dty0(end))
                dy
                idn=find(~isnan(dty));
                lenv=length(idn)-1;
                dtvv1=dtv1(idn,dy); dtvv2=dtv2(idn,dy); dtvv3=dtv3(idn,dy);
                dtvv4=dtv4(idn,dy); dtvv5=dtv5(idn,dy); dtvv6=dtv6(idn,dy);
                dtv1(idn,dy)=NaN; dtv2(idn,dy)=NaN; dtv3(idn,dy)=NaN; 
                dtv4(idn,dy)=NaN; dtv5(idn,dy)=NaN; dtv6(idn,dy)=NaN;

                dtv1(end-lenv:end,dy)=dtvv1; dtv2(end-lenv:end,dy)=dtvv2; dtv3(end-lenv:end,dy)=dtvv3;
                dtv4(end-lenv:end,dy)=dtvv4; dtv5(end-lenv:end,dy)=dtvv5; dtv6(end-lenv:end,dy)=dtvv6;              
            end
        end
        dt3(:,1)=reshape(dtv1,ldt,1); dt3(:,2)=reshape(dtv2,ldt,1); dt3(:,3)=reshape(dtv3,ldt,1);
        dt3(:,4)=reshape(dtv4,ldt,1); dt3(:,5)=reshape(dtv5,ldt,1); dt3(:,6)=reshape(dtv6,ldt,1);
        dtvv=zeros(size(dtv1)); dtvv(~isnan(dtv1))=1;
        dt2=dt3;

    catch ME
        display(ME);
        display(ME.stack);
        display('Error occured while reading in file.');
        continue;
    end

    % Check if the file is empty
    if isempty(dt2)
        display('File is empty');
        continue;
    end
    
    output_path = strcat(read_dir, file_name, '.mat');
    display(replace('Saving %s', '%s', output_path));
    save(output_path,'dt2','fs');
end
display('COMPLETE');
exit(0);
