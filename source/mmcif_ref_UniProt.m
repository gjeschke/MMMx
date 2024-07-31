function reference = mmcif_ref_UniProt(fname)

reference = [];

fid = fopen(fname,'rt');
if isempty(fid) || fid == -1
    return
end
looping = false;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if length(tline) >= length('loop_') && strcmpi(tline(1:length('loop_')),'loop_')
        looping = true;
    end
    if length(tline) >= length('#') && strcmpi(tline(1:length('#')),'#')
        looping = false;
    end
    if length(tline) >= length('_struct_ref_seq.align_id')...
            && strcmpi(tline(1:length('_struct_ref_seq.align_id')),'_struct_ref_seq.align_id')
        if looping
            while tline(1) == '_'
                tline = fgetl(fid);
            end
            while tline(1) ~= '#'
                args = split(tline);
                chain = args{4};
                if double(chain(1)) < double('A') % guard against wrong chain IDs in mmcif
                    chain = char(double(chain(1))+16);
                    if double(chain) < 65 % guard against chain number zero 
                        chain = 'Z';
                    end
                end
                reference.(chain).UniProt = args{9};
                reference.(chain).seq_align_begin = str2double(args{5});
                reference.(chain).seq_align_end = str2double(args{7});
                reference.(chain).db_align_begin = str2double(args{10});
                reference.(chain).db_align_end = str2double(args{12});
                tline = fgetl(fid);
            end
            looping = false;
        else
            for l = 1:3
                tline = fgetl(fid);
            end
            args = split(tline);
            chain = args{2};
            if double(chain(1)) < double('A') % guard against wrong chain IDs in mmcif
                chain = char(double(chain(1))+16);
            end
            if double(chain) < 65 % guard against chain number zero
                chain = 'Z';
            end
            tline = fgetl(fid);
            args = split(tline);
            reference.(chain).seq_align_begin = str2double(args{2});
            fgetl(fid);
            tline = fgetl(fid);
            args = split(tline);
            reference.(chain).seq_align_end = str2double(args{2});
            fgetl(fid);
            tline = fgetl(fid);
            args = split(tline);
            reference.(chain).UniProt = args{2};
            tline = fgetl(fid);
            args = split(tline);
            reference.(chain).db_align_begin = str2double(args{2});
            fgetl(fid);
            tline = fgetl(fid);
            args = split(tline);
            reference.(chain).db_align_end = str2double(args{2});
        end
    end
end
fclose(fid);