classdef data_1P < handle

    properties
        raw_data
        meta
    end

    properties (Dependent)
        
    end

    properties (Access = private)
        rfpCache
        rfp_HD_Cache
        gfpCache
        gfp_HD_Cache
        HbR_Cache
        HbO_Cache
        HbT_Cache
    end

    methods
        function obj = data_1P(path)
            % initialize properties

            data_path = fullfile(path, 'data.bin');
            meta_path = fullfile(path, 'meta.json');
            
            if exist(data_path, 'file') && exist(meta_path, 'file')
                % load metadata from .json
                fid = fopen(meta_path, 'r');
                meta = fread(fid, inf, '*char');
                fclose(fid);
                obj.meta = jsondecode(meta');
                
                % load raw 1P data from .bin
                fid = fopen(data_path, 'r');
                data = fread(fid, inf, obj.meta.dtype);
                fclose(fid);
    
                if obj.meta.order == 'C'
                    order = obj.meta.shape(end:-1:1);
                else
                    order = obj.meta.shape;
                end
    
                obj.raw_data = reshape(data,order');
                obj.raw_data = permute(obj.raw_data, [2, 1, 4, 3]);
            end

            obj.raw_data = rot90(obj.raw_data, obj.meta.rotation / 90);
        end
        
        % calculate dF / F
        function val = rfp(obj)
            if isempty(obj.rfpCache)
                channel_idx = find(obj.meta.channel_order == 565);
                obj.rfpCache = obj.raw_data(:,:,:,channel_idx) ./ mean(obj.raw_data(:,:,:,channel_idx), 3) - 1;
            end
            val = obj.rfpCache;
        end

        function val = gfp(obj)
            if isempty(obj.gfpCache)
                channel_idx = find(obj.meta.channel_order == 470);
                obj.gfpCache = obj.raw_data(:,:,:,channel_idx) ./ mean(obj.raw_data(:,:,:,channel_idx), 3) - 1;
            end
            val = obj.gfpCache;
        end

        % estimate hemodynamics
        function val = HbO(obj)
            if isempty(obj.HbO_Cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);
    
                [obj.HbR_Cache, obj.HbO_Cache] = f_calcHb(obj.raw_data(:,:,:,channel_HD1), obj.raw_data(:,:,:,channel_HD2));
            end
            val = obj.HbO_Cache;
        end

        function val = HbR(obj)
            if isempty(obj.HbR_Cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);
    
                [obj.HbR_Cache, obj.HbO_Cache] = f_calcHb(obj.raw_data(:,:,:,channel_HD1), obj.raw_data(:,:,:,channel_HD2));
            end
            val = obj.HbR_Cache;
        end

        function val = HbT(obj)
            if isempty(obj.HbT_Cache)
                obj.HbT_Cache = obj.HbR + obj.HbO;
            end
            val = obj.HbT_Cache;
        end

        % apply hemodynamic correction
        function val = rfp_HD(obj)
            if isempty(obj.rfp_HD_Cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);

                tmpHbRed = obj.raw_data(:,:,:,channel_HD2) ./ mean(obj.raw_data(:,:,:,channel_HD2),3);
                tmpHbGreen = obj.raw_data(:,:,:,channel_HD1) ./ mean(obj.raw_data(:,:,:,channel_HD1),3);
                obj.rfp_HD_Cache = (obj.rfp + 1) ./ (tmpHbRed.^0.8.*tmpHbGreen.^0.4) - 1;
            end
            val = obj.rfp_HD_Cache;
        end

        function val = gfp_HD(obj)
            if isempty(obj.gfp_HD_Cache)
                tmpWL = [470 515];
                tmpExtinction = f_GetExtinctions(tmpWL);     % in cm
                tmpPathEx = f_pathlengths(tmpWL(1),0.4)/2;       % pathlengths returns in cm
                tmpPathEm = f_pathlengths(tmpWL(2),0.4)/2;
                tmpMuaEx = (tmpExtinction(1,1).*obj.HbO) + (tmpExtinction(1,2).*obj.HbR);
                tmpMuaEm = (tmpExtinction(2,1).*obj.HbO) + (tmpExtinction(2,2).*obj.HbR);
                
                obj.gfp_HD_Cache = (obj.gfp + 1)./exp(-(tmpMuaEx.*tmpPathEx + tmpMuaEm.*tmpPathEm))-1;
            end
            val = obj.gfp_HD_Cache;
        end

        % processing functions

    end

end