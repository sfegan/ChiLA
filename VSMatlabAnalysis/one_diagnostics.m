%
% Crunch one H5 files to diagnostics
%
function one_diagnostics(varargin)
  
  if(nargin<2)
    fprintf(2,'one_diagnostics: need a filename!\n');
    return
  end
  name=varargin{1};
  diagnostics_dat=varargin{2};
  fprintf(1,'Loading %s...\n',name);
  try
    data=loadh5(name);
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
