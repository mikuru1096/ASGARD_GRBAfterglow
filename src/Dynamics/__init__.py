from src.Dynamics.Dynamics_reverse import dynamics_reverse
from src.Dynamics.Dynamics_forward import dynamics_forward


__all__ = [
           "dynamics_reverse",
           "dynamics_forward"
          ]
          
__doc__ = '''
------------------------------------------------------    
|        
| 包含正向和反向激波同时存在,及只考虑正向激波存在的两种情况.  
|                                                    
| dynamics_forward                                   
| 正向激波:
|            
|  求解huang/zhang/peer三种不同的动力学模型,
|  求解算法为自适应步长四阶龙格库塔方法.
|            
| dynamics_reverse
| 反向激波:
|              
|  求解yan的动力学方法,
|  该方法对应纯火球,不考虑喷流的磁化因子.
|
|  当反向激波穿越完毕时,动力学解回到只存在正向激波的情况.
|
------------------------------------------------------
          '''
