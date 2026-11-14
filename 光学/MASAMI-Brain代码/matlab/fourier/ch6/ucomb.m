function[out]=ucomb(x)
% ===== 附录 C.5:ucomb 单位样本"梳"函数(逐行忠实 + 中文注释) =====
% 书附录 C.5 列出的函数(原书第 1-12 行):在 x 为整数值的位置放 1。
% 先把输入四舍五入到 10^-6(第六位小数),再用 rem(x,1)==0 判整数,
% 这样小的舍入误差不会影响 delta 的落位(书 6.4.2 节原话)。
% 使用前提:坐标必须让每个需要的 delta 位置上都有样本——
%   例如 grating_sqr 里 x1/P 的步距 P/dx1 必须是整数(习题 6.6)。
% 注意:输出是 logical(0/1),乘法时 MATLAB 会自动当 0/1 双精度用;
%   Python 版在 python/fourier/ch6/ch6_transmittance.py 的 ucomb()。
% ------------- 以下为书原代码(附录 C.5 第 1-12 行) -------------
% unit sample "comb" function
%
% sequence of unit values for x=integer value
% round is used to truncate roundoff error
%
x=round(x*10^6)/10^6;   %round to 10^6ths place
out=rem(x,1)==0;        %place 1 in out where rem = 0
end
