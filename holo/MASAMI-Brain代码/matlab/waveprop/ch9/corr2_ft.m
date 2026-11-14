function [mcf rgrid] = corr2_ft(U1, U2, mask, delta)
% 复场互相干函数的 FFT 估计（书 Listing 3.10 的复现；ch3 分册 PDF 缺失，
% 按 9.5.5 节/Listing 9.8 的用法重写，与 Python 版数值一致）
%
%   MCF(rho) = sum_r U1(r) U2*(r + rho) mask(r) mask(r + rho)
%
% N×N 循环相关实现：口径 mask 直径（D2）远小于网格边长（N*delta）时，
% 不存在"两点同在口径内又相距近一屏"的卷绕对，循环相关即线性相关。
% 返回 fftshift 后的居中数组（零间隔在 (N/2+1, N/2+1)，与书一致）。
% 用法（Listing 9.8）：MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, deltan);
% 第二个输出 rgrid 是各像素对应的间隔模长 [m]（供径向分箱平均）。

A = fft2(U1 .* mask);
B = fft2(U2 .* mask);             % 共轭在频域做：相关定理 IFT[conj(A).*B]
mcf_full = ifft2(conj(A) .* B);   % wrap 序：C(m) = sum_j U1(j)U2*(j+m)
mcf = fftshift(mcf_full);         % 零间隔居中
N = size(U1, 1);
[ii jj] = ndgrid((0:N-1) - N/2, (0:N-1) - N/2);
rgrid = delta * sqrt(ii.^2 + jj.^2);
end
