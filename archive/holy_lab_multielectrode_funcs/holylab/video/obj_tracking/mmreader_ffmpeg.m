classdef mmreader_ffmpeg < handle
% mmreader_ffmpeg: similar to mmreader in matlab
% See also: mmreader.
   properties(SetAccess = private)
        Name            % Name of the file to be read.
        Path            % Path of the file to be read.
        Duration        % Total length of file in seconds.
   end
   
   properties(GetAccess='public', SetAccess='private', Dependent)
        CurrentTime     % The time of the current frame
        CurrentFrameNumber    % The frame number of the current frame
        NumberOfFrames  % Total number of frames in the video stream. 
        VideoFormat     % Video format as it is represented in MATLAB.
        BitsPerPixel    % Bits per pixel of the video data.
        Height          % Height of the video frame in pixels.
        Width           % Width of the video frame in pixels.
   end
   
   properties(Dependent)
        FrameRate       % Frame rate of the video in frames per second.
   end
   
   properties(Access='private', Hidden)
        fid % ffmpeg id
        filename
        nframes
   end

   methods(Access='public')
      function obj = mmreader_ffmpeg(fileName,varargin)
        % Syntax:
        %   mov = mmreader_ffmpeg('filename','property1','value1',...)
        % This creates a movie reader object. See the various properties
        % for options (and see help for mmreader).
        % One special option: 'calibrate', which if the property is true,
        % will calibrate the frame rate of the movie.
        obj.fid=ffmpeg('open', fileName);
        if(obj.fid==-1)
          error('mmreader_ffmpeg: failed to open file');
        end
        obj.filename=fileName;
        obj.Name=replace_parent_dir(fileName, ''); % TODO
        obj.Path=replace_filename(fileName, '');   % TODO
        obj.nframes=-1;
        % obj.FrameRate = ffmpeg('getFrameRate', obj.fid);
        if ~isempty(varargin)
          for i = 1:2:length(varargin)
            if strcmpi('calibrate',varargin{i})
              if (varargin{i+1})
                fps = calibrateFrameRate(obj);
                ffmpeg('forceFrameRate', obj.fid, fps);
              end
            else
              obj.(varargin{i}) = varargin{i+1};
            end
          end
        end
      end
      
      function fps=calibrateFrameRate(obj)
         fps=-1;
         oldPos=ffmpeg('getCurTime',obj.fid);
         if(oldPos<0)
            ffmpeg('decode', obj.fid);
            oldPos=ffmpeg('getCurTime',obj.fid);
         end
         nMaxTries=3; % NOTE: hard-coded
         nFrames=10; % NOTE: hard-coded
         times=zeros(1, nFrames+1);
         suc=false;
         for idxTry=1:nMaxTries
            for idxFrame=1:nFrames+1
               ffmpeg('decode', obj.fid);
               times(idxFrame)=ffmpeg('getCurTime', obj.fid);
            end
            intervals=diff(times);
            spf=mean(intervals);
            spf_std=std(intervals);
            if(spf_std/spf<0.001) % NOTE: 0.1%. Hard-coded
               suc=true;
               fps=1/spf;
               break;
            end
         end
         % ok, now try to go back to the old position
         obj.seekTime(oldPos);
%          ffmpeg('seekPerfectFrame', obj.fid, oldPos);
%          if(ffmpeg('getCurTime',obj.fid)<0)
%             ffmpeg('decode', obj.fid);
%          end
         if(~suc)
            error('mmreader/calibrateFrameRate: failed the calibration'); % NOTE: should reflect it in the rtn value?
         end
      end
   end

   methods
        function value=get.FrameRate(obj)
            value= ffmpeg('getFrameRate', obj.fid);
        end
        function set.FrameRate(obj, value)
            ffmpeg('forceFrameRate', obj.fid, value);
        end
        function value = get.Duration(obj)
            value = ffmpeg('getDuration', obj.fid);
        end
        function value = get.BitsPerPixel(obj)
            value = 24;
        end
        function value = get.Height(obj)
            value = ffmpeg('getHeight', obj.fid);
        end
        function value = get.NumberOfFrames(obj)
            updateNframes(obj);
            value=obj.nframes;
        end
        function value = get.VideoFormat(obj)
            value = 'RGB24';
        end
        function value = get.Width(obj)
            value = ffmpeg('getWidth', obj.fid);
        end
        function t = get.CurrentTime(obj)
          t = ffmpeg('getCurTime', obj.fid);
        end
        function n = get.CurrentFrameNumber(obj)
          n = round(ffmpeg('getCurTime',obj.fid)*obj.FrameRate);
        end
        
        function result=read(obj, varargin)
           updateNframes(obj);
           if(nargin==1)
              range=[1 obj.nframes];
           elseif(nargin>2)
              error('mmreader_ffmpeg/read: please specify 0 or 1 arg');
           else
              range=varargin{1};
           end
           if(isscalar(range))
              if(isinf(range))
                 range=obj.nframes;
              end
              range=[range range];
           end
           if(length(range)~=2)
              error('mmreader_ffmpeg/read: the frame range must be a scalar or 1x2 vector');
           end
           if(isinf(range(2)))
              range(2)=obj.nframes;
           end
           seek(obj, range(1));
           result=cell(1, diff(range)+1);
           for idx=range(1):range(2)
              result{idx-range(1)+1}=ffmpeg('getFrame', obj.fid);
              if(idx~=range(2))
                 % don't move after read the last frame
                 ffmpeg('decode', obj.fid);
              end
           end % for,
           
           if(length(result)==1)
              result=result{1};
           else
              result=cat(4, result{:});
           end
        end % read()
        function result=readNextFrame(obj)
          ffmpeg('decode', obj.fid);
          result = ffmpeg('getFrame', obj.fid);
        end
        function result = readAtTime(obj,t)
          if isscalar(t)
            dt = t - obj.CurrentTime;
            if (abs(dt) < 0.5/obj.FrameRate)
              % Do nothing
            elseif (dt > 0 && dt < 1.5/obj.FrameRate)
              ffmpeg('decode', obj.fid);
            else
              obj.seekTime(t);
            end
            result = ffmpeg('getFrame', obj.fid);
          else
            result = cell(1,length(t));
            for frameIndex = 1:length(t)
              obj.seekTime(t);
              result{frameIndex} = ffmpeg('getFrame', obj.fid);
            end
          end
        end
   end % methods
    
    methods (Access='public', Hidden)
        function delete(obj)
            ffmpeg('close', obj.fid);
        end
    end
    
    methods (Access='private', Hidden)
       function updateNframes(obj)
            if(obj.nframes>0) 
               return;
            end
            duration=ffmpeg('getDuration', obj.fid);
            obj.nframes = floor(duration*obj.FrameRate-1);
       end
       
       function seek(obj, frameidx)
          curframeidx = obj.CurrentFrameNumber;
          if(curframeidx==frameidx)
             return;
          end
          if (frameidx < curframeidx || frameidx > curframeidx+obj.FrameRate)
            t = frameidx / obj.FrameRate;
            if (t < 0)
              t = 0;
            end
%             ffmpeg('seekPerfectFrame', obj.fid, t);
            ffmpeg('seekFrame', obj.fid, t);
            ffmpeg('decode', obj.fid);
          end
          % Progress incrementally to the desired frame
          while (obj.CurrentFrameNumber < frameidx)
            ffmpeg('decode', obj.fid);
          end
          if (obj.CurrentFrameNumber ~= frameidx)
           warning('mmreader:ffmpeg','Requested frame number %d differs from actual %d',frameidx,obj.CurrentFrameNumber);
         end
       end
       function seekTime(obj, t)
         if (t < 0 || t > obj.Duration)
           error('Cannot seek beyond the edges of the movie');
         end
         tol = 1/obj.FrameRate/2;
         tcur = obj.CurrentTime;
         if (abs(t-tcur) < tol)
           return
         end
         if (t < tcur || t-tcur > 1)
           ffmpeg('seekPerfectFrame', obj.fid, t);
%            ffmpeg('seekFrame', obj.fid, t);
         end
         while (t - obj.CurrentTime > tol)
           ffmpeg('decode', obj.fid);
         end
         if (abs(t-obj.CurrentTime) > tol)
           warning('mmreader:ffmpeg','Requested time %g differs from actual time %g',t,obj.CurrentTime);
         end
       end
    end
end % classdef
   