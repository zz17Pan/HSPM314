classdef ResourceAllocation
    methods (Static)
        function [comm_subarrays, sens_subarrays] = allocate_subarrays()
            try
                cfg = Config;
                sens_subarrays = struct('tx', 1, 'rx', 1);
                comm_subarrays = struct('tx', 2:cfg.Kt, 'rx', 2:cfg.Kr);
                fprintf('���������:\n');
                fprintf('��֪����: tx=%d, rx=%d\n', sens_subarrays.tx, sens_subarrays.rx);
                fprintf('ͨ������: tx=%d-%d, rx=%d-%d\n', comm_subarrays.tx(1), comm_subarrays.tx(end), comm_subarrays.rx(1), comm_subarrays.rx(end));
            catch ME
                fprintf('�������ʧ�ܣ�\n');
                fprintf('������Ϣ��%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
