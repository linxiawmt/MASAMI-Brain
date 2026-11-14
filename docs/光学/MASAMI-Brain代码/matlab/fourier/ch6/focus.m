function[uout]=focus(uin,L,lambda,zf)
% ===== 书 6.2 节:focus 聚焦/发散透过率函数(逐行忠实 + 中文注释) =====
% 书 6.2 节列出的函数(原书第 1-18 行),实现式 (6.9):
%   t_A(x,y) = exp[-j (k/(2*zf)) (x^2 + y^2)]
% zf > 0 会聚(焦点在源面前方 z = zf 处);zf < 0 发散(虚焦点)。
% 这是与菲涅耳脉冲响应 h(式 4.23)同形的相位 chirp,只是指数符号为负
% ——h 让波前发散开去,focus 把波前向焦点收拢。
% 采样限制(式 6.11):|zf| >= (D1/2/lambda) * (1/(2*dx) - B1)^(-1),
% 即聚焦 chirp 会加宽带宽,大约增加 D1/(2*lambda*|zf|)(式 6.10)。
% MATLAB/Python 差异:Python 版在 python/fourier/ch6/ch6_transmittance.py
%   的 focus();X.^2 在 numpy 里写作 X**2。
% ------------- 以下为书原代码(第 1-18 行) -------------
% converging or diverging phase-front
% uniform sampling assumed
% uin - input field
% L - side length
% lambda - wavelength
% zf - focal distance (+ converge, - diverge)
% uout - output field

[M,N]=size(uin);        %get input field array size
dx=L/M;                 %sample interval
k=2*pi/lambda;          %wavenumber
%
x=-L/2:dx:L/2-dx;       %coords
[X,Y]=meshgrid(x,x);

uout=uin.*exp(-j*k/(2*zf)*(X.^2+Y.^2)); %apply focus
end
