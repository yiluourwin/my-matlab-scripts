v=[7.8 3.6 3.0];
b=300;

filelist=dir;

for i=3:length(filelist)

if isempty(strfind(filelist(i).name,'.SAC'))
continue;
end

[sachd,sacdata]=rsac(filelist(i).name);

t=(0:(length(sacdata)-1)).*sachd.delta;
t=t';

t1=sachd.dist./v(1);
twindow=sachd.dist./v(2:3)-t1+b;

sacdatax=sacdata(t<1000);

tx=t(t<1000);

%twindow=twindow-tx(1);

Fs=1/sachd.delta;
[a,f]=spectra_cal(Fs,sacdatax,[0.01,10],0.05,twindow);

out=cat(2,f,a);
out=out';

fid=fopen([filelist(i).name,'.spec'],'w');

fprintf(fid,'%f %f\n',out);
fclose(fid);

%
twindow=sachd.dist./v(2:3)-sachd.dist./v(3)+b;

[a,f]=spectra_cal(Fs,sacdatax,[0.01,10],0.05,twindow);

out=cat(2,f,a);
out=out';

fid=fopen([filelist(i).name,'.noise'],'w');

fprintf(fid,'%f %f\n',out);
fclose(fid);

end
