whereami = pwd;
response = cx_command(sprintf('cd %s',whereami),options);
% response = ChimeraX('cd g:\MMM_develop\ChimeraInterface');
% fprintf(1,'CX>%s\n',response);
% [id,entity] = pdbload('2LZM');
response = ChimeraX('open 2lzm');
fprintf(1,'CX>%s\n',response);
response = ChimeraX(sprintf('open %s name loc131R1',cubename));
fprintf(1,'CX>%s\n',response);
response = ChimeraX('color #2 crimson');
fprintf(1,'CX>%s\n',response);
response = ChimeraX('transparency #2 50');
fprintf(1,'CX>%s\n',response);

