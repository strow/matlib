function varargout = getdata_opendap(url, field)
% Usage: vars = getdata_opendap(url, fields)
% A low level function call or accessing opendap data sets.
%
% url - the url of the opendap data set to read from
% field - field name and dimensions of the data to retrieve
%
%  For example:
% [r s] = getdata_opendap(['http://airscal2u.ecs.nasa.gov/opendap/Aqua_AIRS_Level1/' ...
%  'AIRIBRAD.005/2012/177/AIRS.2012.06.25.001.L1B.AIRS_Rad.v5.0.21.0.G13013010558.hdf'], ...
%  'radiances[0:1:134][44:1:45][0:1:2377],scanang[0:1:134][44:1:45]');
%


if(length(url) == 1)
	url = url{1};
end

t = mktemp();
disp([url '.dods?' field])
urlwrite([url '.dods?' field],t);

% Preparse out the field dimensions to check validity
f_struct = regexp([',' field],'[,](?<var>\w+)(?<span>\[[^,]+)','names');
for i = 1:length(f_struct)
    in_names{i} = f_struct(i).var;
    in_dims{i} = cellfun(@(v) length(eval(v{1})),regexp(f_struct(i).span,'\[([0-9:]+)\]','tokens'));
end

var_dims = {};
var_names = {};
var_type = {};
var_store = {};
count = 1;
h=fopen(t,'r','ieee-be');
for i=1:100
    str = fgets(h);
    if(strcmp(['Data:'],str(1:5)))
        break
    end
    
    % clear out spaces
    s = str(str ~= ' ');
    if ~ischar(s)
        continue
    end

    struct_re = regexp(str,' (\w*) (\w*)__0;','tokens');
    if length(struct_re) > 0
        var_type{count} = struct_re{1}{1};
        var_names{count} = struct_re{1}{2};
        var_store{count} = 'struct';
        for j = 1:length(in_names)
            if strcmp(var_names{count},in_names{j})
                var_dims{count} = in_dims{j};
            end
        end
        count = count + 1;
    end

    name_re = regexp(str,'^ *(\w*) (\w*)[','tokens');
    if length(name_re) > 0
        var_type{count} = name_re{1}{1};
        var_names{count} = name_re{1}{2};
        var_store{count} = 'block';

        dims = cellfun(@(i) str2num(i{1}),regexp(str,'= (\d)*]','tokens'));
				if count <= length(in_names)
          for j = 1:length(in_names)
            if strcmp(var_names{count},in_names{j})
                if ~isequal(dims, in_dims{j})
                    error(['dimension mismatch ' num2str(var_dims{count}) ' and ' num2str(dims)])
                end
                var_dims{count} = dims;
            end
          end
        else
					var_dims{count} = dims;
        end

        count = count + 1;
    end

end

for i = 1:length(var_store)

    %if i > 1
    %    fseek(h,-24,0);
    %    reshape(fread(h,48,'uint8'),8,[])
    %    fseek(h,-24,0);
    %end

    if strcmp(var_store{i},'block')
        dims = fread(h,2,'int32');
        if dims(1) ~= dims(2)
          % A weird thing in OpenDAP, that two 0 bytes get inserted in the binary stream...
          %   to avoid running into trouble let's skip these.  This normally shows up after type 10.
          %disp('padding')
          fseek(h,-6,0);
          dims = fread(h,2,'int32');
        end
        if dims(1) ~= dims(2)
            error(['getdata_opendap:  Dimensions declared in binary header in OpenDAP field ''' var_names{i} ''' is not consistent.']);
        end

        if prod(var_dims{i}) ~= dims(1)
            error(['getdata_opendap:  Dimensions on OpenDAP request field ''' var_names{i} ''' is not consistent with what was returned.']);
        end
    end

    if length(var_dims{i}) == 1; var_dims{i} = [var_dims{i} 1]; end %make 2d

    % Parse out the variable types in_cast = fread format, typesize = binning size in struct data, type = cast format in struct data
    if(strcmpi(var_type{i},'Float64'))
        in_class = 'float64';  typesize = 8;  type = 'double';
    elseif(strcmpi(var_type{i},'Float32'))
        in_class = 'float32';  typesize = 4;  type = 'single';
    elseif(strcmpi(var_type{i},'Byte'))
        in_class = 'uint8';    typesize = 4;  type = 'uint32';
        % normally two bytes are padded here (very odd!)
    elseif(strcmpi(var_type{i},'Int8'))
        in_class = 'int8';     typesize = 4;  type = 'int32';
    elseif(strcmpi(var_type{i},'Int16'))
        in_class = 'int32';    typesize = 4;  type = 'int32';
    elseif(strcmpi(var_type{i},'Int32'))
        in_class = 'int32';    typesize = 4;  type = 'int32';
    elseif(strcmpi(var_type{i},'Int64'))
        in_class = 'int64';    typesize = 8;  type = 'int64';
    end

    if strcmp(var_store{i},'block')
        varargout{i} = reshape(fread(h,dims(1),in_class),var_dims{i});
    elseif strcmp(var_store{i},'struct')
        if ~isequal(fread(h,4,'uint8')',[90 0 0 0])
            fseek(h,-2,0);
            if ~isequal(fread(h,4,'uint8')',[90 0 0 0])
                 error(['getdata_opendap:  Field ''' var_names{i} ''' binary start is not a structure format.']);
            end
        end
        fseek(h,-4,0);

        data = reshape(fread(h,var_dims{i}(1)*(4+typesize),'uint8'),[],var_dims{i}(1));
        %reshape(data(4+(1:4),:),1,[])
        if ~isequal(unique(data((1:4),:)','rows'),[90 0 0 0])
            error(['getdata_opendap:  Field ''' var_names{i} ''' binary is not a structure.']);
        end
        varargout{i} = swapbytes(typecast(reshape(uint8(data(4+(1:typesize),:)),1,[]),type));
        fread(h,4,'uint8');
        %keyboard
        %varargout(i) = {reshape(fread(h,var_dims{i}(1)*8,'uint8'),var_dims{i},[])};
    end
    %if ndims(varargout{i}) < 3
    %  plot(varargout{i});title(var_names{i})
    %  pause
    %end

end
%dat = fread(h,inf);
fclose(h);
delete(t);

if i == length(in_names)
  % Order the data in the request order for parsing
  [a b] = ismember(in_names,var_names);
  varargout = varargout(b);
end
