% function to read initial state from file
function [etp_c]=read_etp_new(etp_file)
fid=fopen(etp_file,'r');
comment=fgetl(fid); %just for reading the right lines in the subsequent section
f_con=str2num(fgetl(fid)); %reads conversion factors
etp_c=[];
while feof(fid)==0 %loop as long as values are in the etp_file
help_array=str2num(fgetl(fid)); %reads the next line of the etp_file and stores the values for time and etp temporarily in a kind of 'container'
help_array(1)=f_con(1)*help_array(1); %convert units to si units in this container for 1) time
help_array(2)=f_con(2)*help_array(2); %convert units to si units in this container for 2) etp
etp_c=[etp_c;help_array]; %stores the values of the container in the final etp_c matrix
end
fclose(fid);
