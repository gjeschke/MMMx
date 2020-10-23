entity0 = get_pdb('2lzm');

whereami = pwd;
response = cx_command(sprintf('cd %s',whereami));
if strcmpi(response,'### FAILED ###')
    fprintf(2,'ChimeraX not available.\n');
else
    fprintf(1,'CX> %s\n',response);
    [entity,exception,id] = cx_get_pdb('1a6m');
    fprintf(1,'CX> %s\n',response);
end
disp('Aber hallo');


