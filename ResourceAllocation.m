classdef ResourceAllocation
    methods (Static)
        function [comm_subarrays, sens_subarrays] = allocate_subarrays()
            try
                cfg = Config;
                sens_subarrays = struct('tx', 1, 'rx', 1);
                comm_subarrays = struct('tx', 2:cfg.Kt, 'rx', 2:cfg.Kr);
                fprintf('子阵分配结果:\n');
                fprintf('感知子阵: tx=%d, rx=%d\n', sens_subarrays.tx, sens_subarrays.rx);
                fprintf('通信子阵: tx=%d-%d, rx=%d-%d\n', comm_subarrays.tx(1), comm_subarrays.tx(end), comm_subarrays.rx(1), comm_subarrays.rx(end));
            catch ME
                fprintf('子阵分配失败：\n');
                fprintf('错误信息：%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
