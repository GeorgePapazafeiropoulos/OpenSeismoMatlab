function  [RecordName_all, dt_all, numPoint_all, PGA_all, Acc_all] =...
fn_ParseAT2File2(GMFileLocation, records)

% Non-scaling
SF_all = ones(1, length(records{1}));


%% Generate AT2 File Data
global recd_time recd_acc       % earthquake record
gravity=9.81;                  % gravity accleration, m/sec^2

for r = 1:length (records{1})
     SF = SF_all(r);
    %	Read earthquake acceleration record from file within Peer 2019 format
    RecordName0 = cellstr(records{1}(r));
    RecordName = strcat(GMFileLocation,RecordName0);
    filename=cellstr(RecordName);
    fid2=fopen(filename{1}, 'r');
    
    %   Skip the first three lines
    header1 = fgetl(fid2);
    header2 = fgetl(fid2);
    header3 = fgetl(fid2);
    
    %    Catch the number of data point (NPTS) and time interval (dt) in row 4
    if( header1(1:4) ~= 'PEER')
        display('The earthquake record is not in PEER format!')
    else
        temp2 = textscan(fid2, '%*s %f %*s %*s %f %*[^\n]', 1);
    end
    
    %	Catch the accel record data in 'data' cell, and transform the amlitude to "recor_acc" var,
    %	transfrom the record time to "recd_time" var.
    numPoint=temp2{1,1};     % number of data points
    dt=temp2{1, 2};      % time step, sec
    fmt='%f %f %f %f %f';   % format of each line in the record file
    numRow = floor( numPoint / 5 );   % number of rows to read
    numVoid = numPoint - 5*numRow;
    data = textscan(fid2, fmt, numRow, 'CollectOutput', 1);
    if numVoid~=0
        fmt2 = '%f';
        for i =1:numVoid-1
            fmt2 = [fmt2, ' %f'];
        end
        data0 = textscan(fid2, fmt2, 1,'CollectOutput', 1);
    end
    fclose(fid2);
    record=data{:};                     %unit: g
    record=record';
    recd_acc=reshape(record, [], 1);    %unit: g, in one-col vector
    if numVoid~=0
        %      recd_acc=[recd_acc;(data0{1}.* Scaling(r))'];
        recd_acc=[recd_acc;data0{1}'];
    end
    recd_acc = recd_acc * SF;	% unit: g
    
   % Writhe the record accleartion into a cell
    
    PGA=max(abs(recd_acc));             %unit: g
    
    recd_acc=gravity.*recd_acc;         %unit: m/s^2
    Acc_all{r} = recd_acc;   
 
    % Write Input Info. into a Matrix
    RecordName_all(r) = RecordName0;
    dt_all(r) = dt;
    numPoint_all(r) = numPoint;
    PGA_all(r) =PGA ;

end

%
%__________________________________________________________________________
% Copyright (c) 2022
%     Chenhao Wu
%     Southeast University (China)
%     Email: wch@seu.edu.cn
% _________________________________________________________________________
