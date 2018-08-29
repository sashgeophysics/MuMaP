clear;

cd('~/codes/ROPT/input');
data=load('Hawaii_RF_Results.dat');
n=size(data,1);
ncols=size(data,2);
%# create format string for fprintf. Use repmat to replicate the %8.2f's
fmtString = [repmat('%16.5f',1,ncols),'\n'];


fid=fopen('Hawaii_RF_Results.csv');
n=size(data,1);
for ii=1:n
   fprintf(fid,fmtString,data(ii,:)) ;
end
fclose(fid);