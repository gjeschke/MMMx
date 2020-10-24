whichone = 'color';
value = 'darkgreen';

clear options
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% id = cx_get_pdb('1a6m',options);
% id = cx_get_pdb('2ad9',options);
% id = cx_get_pdb('2lzm',options);

% address = '(A)57';
address = '(A){3}';
% address = '(A)154.FE';
% address = '(A)156.S';
attribute = cx_set_atom(id,address,whichone,value);
