function G = ft2(g, delta)
% 书 Schmidt PM199 Listing 2.5：定标正确的二维 DFT（复现，ch9 依赖）
% G(fx, fy) = ∫∫ g(x,y) e^{-i2π(fx x + fy y)} dx dy 的离散近似：
% 双侧 fftshift 把零频移到中心，乘 delta^2 对应 dx dy 的面积元。
% MATLAB/Python 差异：numpy 侧同款实现见 fourier/common/base.py 的
% ft2_centered；两侧数值行为一致（居中排列，零频在 N/2+1 / N//2）。
G = fftshift(fft2(fftshift(g))) * delta^2;
end
