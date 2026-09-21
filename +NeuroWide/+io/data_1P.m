classdef data_1P < handle

    properties
        raw_data
        meta
    end

    properties (Dependent, Hidden)
        rfp
        rfp_HD
        gfp
        gfp_HD
        HbO
        HbR
        HbT
    end

    properties (Access = private)
        rfp_cache
        rfp_HD_cache
        gfp_cache
        gfp_HD_cache
        HbO_cache
        HbR_cache        
    end

    methods
        function obj = data_1P(path)
            % initialize properties
            data_path = fullfile(path, 'data.bin');
            meta_path = fullfile(path, 'meta.json');
            
            if exist(data_path, 'file') && exist(meta_path, 'file')
                % load metadata from .json
                % fid = fopen(meta_path, 'r');
                % meta = fread(fid, inf, '*char');
                % fclose(fid);
                obj.meta = jsondecode(fileread(meta_path));
                
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
                obj.raw_data = permute(obj.raw_data, [2,1,4,3]);
            end

            obj.raw_data = rot90(obj.raw_data, obj.meta.rotation / 90);
        end

        function disp(obj)
            fprintf(['Data loader class for mesoscale imaging data\n\n' ...
                'Attributes:\n' ...
                'rfp: red fluorescence delta F / F\n' ...
                'gfp: green fluorescence delta F / F\n' ...
                'rfp_HD: hemodynamic corrected red fluorescence\n' ...
                'gfp_HD: hemodynamic corrected green fluorescence\n' ...
                'HbO: estimated oxygenated hemoglobin\n' ...
                'HbR: estimated deoxygenated hemoglobin\n' ...
                'HbT: estimated total hemoglobin\n']);
        end
        
        % calculate dF / F
        function val = get.rfp(obj)
            if isempty(obj.rfp_cache)
                channel_idx = find(obj.meta.channel_order == 565);
            obj.rfp_cache = obj.raw_data(:,:,:,channel_idx) ./ mean(obj.raw_data(:,:,:,channel_idx), 3) - 1;
            end
            val = 100 * obj.rfp_cache;
        end

        function val = get.gfp(obj)
            if isempty(obj.gfp_cache)
                channel_idx = find(obj.meta.channel_order == 470);
                obj.gfp_cache = obj.raw_data(:,:,:,channel_idx) ./ mean(obj.raw_data(:,:,:,channel_idx), 3) - 1;
            end
            val = 100 * obj.gfp_cache;
        end

        % estimate hemodynamics
        function val = get.HbO(obj)
            if isempty(obj.HbO_cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);
    
                [obj.HbR_cache, obj.HbO_cache] = NeuroWide.process.calcHb(obj.raw_data(:,:,:,channel_HD1), obj.raw_data(:,:,:,channel_HD2));
            end
            val = 1e6 * obj.HbO_cache;
        end

        function val = get.HbR(obj)
            if isempty(obj.HbR_cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);
    
                [obj.HbR_cache, obj.HbO_cache] = NeuroWide.process.calcHb(obj.raw_data(:,:,:,channel_HD1), obj.raw_data(:,:,:,channel_HD2));
            end
            val = 1e6 * obj.HbR_cache;
        end

        function val = get.HbT(obj)
            val = obj.HbR + obj.HbO;
        end

        % apply hemodynamic correction
        function val = get.rfp_HD(obj)
            if isempty(obj.rfp_HD_cache)
                channel_HD1 = find(obj.meta.channel_order == 525);
                channel_HD2 = find(obj.meta.channel_order == 625);
    
                tmpHbRed = obj.raw_data(:,:,:,channel_HD2) ./ mean(obj.raw_data(:,:,:,channel_HD2),3);
                tmpHbGreen = obj.raw_data(:,:,:,channel_HD1) ./ mean(obj.raw_data(:,:,:,channel_HD1),3);
                obj.rfp_HD_cache = ((obj.rfp / 100 + 1) ./ (tmpHbRed.^0.8.*tmpHbGreen.^0.4) - 1);
            end
            val = 100 * obj.rfp_HD_cache;
        end

        function val = get.gfp_HD(obj)
            if isempty(obj.gfp_HD_cache)
                tmpWL = [470,515];
                tmpExtinction = NeuroWide.process.getExtinctions(tmpWL);     % in cm
                tmpPathEx = NeuroWide.process.pathlengths(tmpWL(1),0.4)/2;       % pathlengths returns in cm
                tmpPathEm = NeuroWide.process.pathlengths(tmpWL(2),0.4)/2;
                tmpMuaEx = (tmpExtinction(1,1).*obj.HbO / 1e6) + (tmpExtinction(1,2).*obj.HbR / 1e6);
                tmpMuaEm = (tmpExtinction(2,1).*obj.HbO / 1e6) + (tmpExtinction(2,2).*obj.HbR / 1e6);
    
                obj.gfp_HD_cache = ((obj.gfp / 100 + 1)./exp(-(tmpMuaEx.*tmpPathEx + tmpMuaEm.*tmpPathEm)) - 1);
            end
            val = 100 * obj.gfp_HD_cache;
        end

        % processing functions

    end

end