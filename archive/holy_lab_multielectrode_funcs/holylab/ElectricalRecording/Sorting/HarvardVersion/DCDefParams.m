function params = DCDefParams(params)
if (nargin < 1)
  params = struct;
end
if (~isfield(params,'dispsnips'))
        params.dispsnips = 50;
end
if (~isfield(params,'ACTime'))
        params.ACTime = 0.1;
end
if (~isfield(params,'NSpikes'))
        params.NSpikes = 2000;
end
if (~isfield(params,'NNoise'))
        params.NNoise = 5000;
end
if (~isfield(params,'BlockSize'))
        params.BlockSize = 5000;
end
if (~isfield(params,'ClustNumOffset'))
        params.ClustNumOffset = 0;
end
