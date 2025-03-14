classdef Config < handle
    properties (Constant)
        c = 3e8                
        f = 300e9             
        lambda = 0.001        
        Nx = 4                
        Nz = 4                
        Kt = 4                
        Kr = 4                
        d = 0.0005           
        d_sub = 0.016        
        k_f = 0.1            
        Beta = 1e7           
        Ts = 0.01            
        SNR = 20             
        fmcw_bandwidth = 1e9;    
        sweep_time = 1e-3;       
        sample_rate = 2e9;       
    end
end
