function [] = convert_svg_pdf(filename)
inkscape='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
svg = [ strrep(pwd, ' ', '\ ') '/' filename ];
pdf = [ svg(1:end-3) 'pdf' ];
command = [inkscape ' ' svg ' --export-pdf=' pdf ];
system(command);
end