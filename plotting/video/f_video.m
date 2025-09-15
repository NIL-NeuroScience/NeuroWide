classdef f_video
    properties
        v
    end

    methods
        function obj = f_video(name,varargin)
            
            p = inputParser;
            addParameter(p,'fs',1);
            parse(p,varargin{:});
            
            obj.v = VideoWriter(name,'Motion JPEG AVI');
            obj.v.FrameRate = p.Results.fs;
            open(obj.v);
        end
        
        function write(obj)
            frame = getframe(gcf);
            writeVideo(obj.v,frame);
        end
        
        function close(obj)
            close(obj.v);
            clear obj;
        end
    end
end