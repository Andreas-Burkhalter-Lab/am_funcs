function varargout=ffmpeg(action, varargin)
% NOTE: the proper order of calling is:
%       open [seek decodeIntoFrame getFrame]... close
% NOTE: after calling decode, the pos is also "seeked" to the next decodable frame

   switch action
      case 'open'
         varargout{1}=mat_ffmpeg(1, varargin{1});
      case 'close'
         mat_ffmpeg(2, varargin{1});
      case 'getCurTime'
         varargout{1}=mat_ffmpeg(3, varargin{1});
      case 'decode' % decode cur frame
         varargout{1}=mat_ffmpeg(4, varargin{1});
      case 'seekFrame'
         varargout{1}=mat_ffmpeg(5, varargin{1}, varargin{2});
      case 'seekPerfectFrame'
         varargout{1}=mat_ffmpeg(6, varargin{1}, varargin{2});
      case 'getHeight'
         varargout{1}=mat_ffmpeg(7, varargin{1});
      case 'getWidth'
         varargout{1}=mat_ffmpeg(8, varargin{1});
      case 'getDuration' % 
         varargout{1}=mat_ffmpeg(9, varargin{1});
      case 'getFrame' % return decoded frame
         data=mat_ffmpeg(10, varargin{1});
         height=mat_ffmpeg(7, varargin{1});
         width =mat_ffmpeg(8, varargin{1});
         data=reshape(data, 3, width, height);
         varargout{1}=permute(data, [3 2 1]);
      case 'getErrorString'
         varargout{1}=mat_ffmpeg(11, varargin{1});
      case 'getFrameRate'
         varargout{1}=mat_ffmpeg(12, varargin{1});
      case 'forceFrameRate'
         mat_ffmpeg(13, varargin{1}, varargin{2});
      otherwise
         error('unknown action');
   end