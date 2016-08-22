function s = readconditionfile(infile)
% function s = readconditionfile(infile[,'2'])

if exist('flag','var'), flag=true; else flag=false; end

fid=fopen(infile,'r');
fgetl(fid);
s1 = 0; s={};
while ~isempty(s1)
    s1 = fscanf(fid,'%d %d %*d %*d %d %d %d',[5 500000]);
    s{end+1} = s1';
end
fclose(fid);
s = cat(1,s{:});

