%
% Crunch all H5 files to diagnostics
%
function all_diagnostics(varargin)
    diagnostics_dat='diagnostics.dat';
    pattern='x*_s2.h5';
    path='';
    if(nargin>0)
        pattern=varargin{1};
        isep=find(pattern=='/',1,'last');
        if(~isempty(isep))
            path=pattern(1:isep);
        end
    end
    files=dir(pattern);
    if(~isempty(files))
      delete(diagnostics_dat);
    end
    for ifile = 1:length(files)
        name=files(ifile).name;
        fprintf(1,'Loading %s%s...\n',path,name);
        try
	    clear data
            data=loadh5(strcat(path,name));
            if(isstruct(data) && isfield(data,'diagnostics'))
                mat_diagnostics(data,[],diagnostics_dat);
            end
        catch
            disp(['Diagnostics for ' name ' failed']);
	    msg=lasterror;
	    disp(msg.message)
	    for istack=1:length(msg.stack)
                disp(sprintf('In %s : %s at %d',msg.stack(istack).file,...
                             msg.stack(istack).name,...
                             msg.stack(istack).line));
            end
        end
    end
end
