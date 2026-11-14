function [xn yn Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)
% ===== 书 Schmidt PM199 Listing 9.1 =====
% 角谱多平面传播（split-step 束传播法，式 9.3）：真空衍射算符 R 与
% 折射（相位屏）算符 T 交替；t 是 N×N×n 的"相互作用"数组（吸收边界
% 或 吸收边界.*exp(i*相位屏)），t(:,:,1) 在源面、t(:,:,idx+1) 在第 idx
% 次部分传播之后。T = 1（全 1 数组）时退化为第 8 章的真空多平面传播。
%
% 原书笔误：式 (9.2) 印成 T = exp(-i*phi)。按时间约定 exp(-i2*pi*nu*t)
% 的抛物线方程 dU/dz = i/(2k)*∇²U + i*k*delta_n*U，折射相位应为
% T = exp(+i*phi)，phi(r) = k*∫delta_n dz；书自己的 Listing 9.8
% （t = sg.*exp(i*phz)）用的也是 + 号。对零均值随机相位屏无统计影响。
%
% 用法：
%   [xn yn Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)
%   Uin:      N×N 源面复场（N 取 2 的幂配合 FFT）
%   wvl:      光波长 [m]
%   delta1:   源面网格间距 [m]
%   deltan:   观测面网格间距 [m]（中间平面按式 8.8 线性插值）
%   z:        部分传播平面位置向量（不含 0），长度 n-1
%   t:        N×N×n 相互作用数组
N = size(Uin, 1);                            % number of grid points
[nx ny] = meshgrid((-N/2 : 1 : N/2 - 1));
k = 2*pi/wvl;                                % optical wavevector
% 超高斯吸收边界（第 8 章式 8.9）：抑制网格边缘卷绕
nsq = nx.^2 + ny.^2;
w = 0.47*N;
sg = exp(-nsq.^8/w^16); clear('nsq', 'w');

z = [0 z];                                   % propagation plane locations
n = length(z);                               % n-1 次部分传播、n 个平面
Delta_z = z(2:n) - z(1:n-1);                 % propagation distances
% 网格间距沿 z 线性插值（式 8.8）：alpha_i = z_i / z_n
alpha = z / z(n);
delta = (1-alpha) * delta1 + alpha * deltan;
m = delta(2:n) ./ delta(1:n-1);              % 相邻平面网格缩放比
x1 = nx * delta(1);
y1 = ny * delta(1);
r1sq = x1.^2 + y1.^2;
Q1 = exp(i*k/2*(1-m(1))/Delta_z(1)*r1sq);    % 源面二次相位（式 8.18 的 Q1）
Uin = Uin .* Q1 .* t(:,:,1);                 % 源面相互作用（相位屏 1）
for idx = 1 : n-1
    % spatial frequencies (of i^th plane)
    deltaf = 1 / (N*delta(idx));
    fX = nx * deltaf;
    fY = ny * deltaf;
    fsq = fX.^2 + fY.^2;
    Z = Delta_z(idx);                        % propagation distance
    % quadratic phase factor（频域啁啾 Q2，与第 6 章 Listing 6.5 同款）
    Q2 = exp(-i*pi^2*2*Z/m(idx)/k*fsq);
    % compute the propagated field：吸收边界 .* 折射 .* [IFT(Q2.*FT(U/m))]
    Uin = sg .* t(:,:,idx+1) ...
        .* ift2(Q2 ...
        .* ft2(Uin / m(idx), delta(idx)), deltaf);
end
% observation-plane coordinates
xn = nx * delta(n);
yn = ny * delta(n);
rnsq = xn.^2 + yn.^2;
Q3 = exp(i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);  % 观测面二次相位
Uout = Q3 .* Uin;
end
