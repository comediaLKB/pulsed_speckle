function [] = remove_diagonal_svg(filename)
fid  = fopen(filename,'r');
f = fread(fid,'*char')';
fclose(fid);
f = regexprep(f,'style="clip-path:url\(#clipPath\d+\);"','');
fid  = fopen(filename,'w');
fprintf(fid,'%s',f);
fclose(fid);
end