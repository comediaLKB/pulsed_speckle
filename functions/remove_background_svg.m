function [] = remove_background_svg(filename)
fid  = fopen(filename,'r');
f = fread(fid,'*char')';
fclose(fid);
f = strrep(f,'fill:white','opacity:0');
fid  = fopen(filename,'w');
fprintf(fid,'%s',f);
fclose(fid);
end