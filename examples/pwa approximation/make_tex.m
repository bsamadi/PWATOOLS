function make_tex(filename)

matlab2tikz([filename '.tikz'],'showInfo', false);
fid = fopen([filename '.tex'],'w');
fprintf(fid,'%s\n','\documentclass{article}','\usepackage{tikz,pgfplots}','\begin{document}','\pagestyle{empty}',['\input{' filename '.tikz}'],'\end{document}');
fclose(fid);
