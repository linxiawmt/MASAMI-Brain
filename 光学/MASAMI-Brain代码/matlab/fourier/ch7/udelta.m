function[out]=udelta(x)
% ===== 附录 C.6:udelta 单位样本 delta(逐行忠实 + 中文注释) =====
% 书附录 C.6 列出的函数(原书第 1-11 行)。
% 在 x=0(六位舍入后)处置 1,用法与 ucomb 同一家族:
%   习题 7.10 里点源 = udelta((X-x0)/dx) .* udelta((Y-y0)/dy),
%   前提是 (x0,y0) 恰好落在采样网格上(holography 例:x0=2 mm、
%   dx=0.02 mm,2/0.02=100 为整数,能精确落位)。
% Python 版在 python/fourier/ch7/ch7_imaging.py 的 udelta()。
% ------------- 以下为书原代码(附录 C.6 第 1-11 行) -------------
%
% unit sample "delta" function
%
% unit value for x=0
% round is used to truncate roundoff error
%
x=round(x*10^6)/10^6;  %round to 10^6ths place
out=x==0;              %place 1 in out where x = 0
end
