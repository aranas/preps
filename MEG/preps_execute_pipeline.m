function preps_execute_pipeline(pipelinename, varargin)

% MOUS_EXECUTE_PIPELINE serves the purpose to make a script executable by qsub.
% supply it with the name of the script that has to be run


if numel(varargin)>0
  for k = 1:numel(varargin)
    if numel(varargin{k})==1
        eval(varargin{k}{1});
    else
        eval([varargin{k}{1},'=varargin{k}{2}']);
    end
  end
end
eval(pipelinename);
