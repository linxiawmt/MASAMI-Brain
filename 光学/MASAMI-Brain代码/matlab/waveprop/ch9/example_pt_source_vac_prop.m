% ===== 书 Schmidt PM199 Listing 9.7 =====
% 例题（9.5.3 节）：真空仿真（先于湍流仿真，两个目的）：
%   1) 不含湍流时结果必须对得上解析预期——ROI 内辐照度均匀、准直后
%      相位平坦（原书图 9.7）；
%   2) 留下真空场 Uvac 供湍流仿真对比（算 Strehl 比、log 振幅扰动
%      chi = ln(|U|/|Uvac|) 都要用它）。
% 前提：先运行 example_pt_source_atmos_setup.m 与
% analysis_pt_source_atmos_samp.m（用到 wvl、k、Dz、D1、R、nscr、d1、d2、N）。
% 依赖：ang_spec_multi_prop.m（Listing 9.1）、sinc（Signal Processing
% Toolbox；无该工具箱可自定义 sinc(x) = sin(pi*x)./(pi*x)）。

delta1 = d1;         % source-plane grid spacing [m]
deltan = d2;         % observation-plane grid spacing [m]
n = nscr;            % number of planes

% coordinates
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1 r1] = cart2pol(x1, y1);

% point source（sinc-Gaussian 点源模型，式 9.3 的 U(r1)）：
%   中央瓣全宽 2*D1 的 sinc 剖面 + 收敛球面波相位 exp(-i*k*r1^2/(2*R))
%   + 高斯包络 exp(-(r1/(4*D1))^2) 抑制旁瓣；/D1^2 定标辐照度量级
pt = exp(-i*k/(2*R) * r1.^2) / D1^2 ...
    .* sinc(x1/D1) .* sinc(y1/D1) ...
    .* exp(-(r1/(4*D1)).^2);

% partial prop planes（n-1 = 10 次部分传播，均匀铺开）
z = (1 : n-1) * Dz / (n-1);

% simulate vacuum propagation（t = 纯吸收边界，无相位屏）
sg = exp(-(x1/(0.47*N*d1)).^16) ...
    .* exp(-(y1/(0.47*N*d1)).^16);
t = repmat(sg, [1 1 n]);
[xn yn Uvac] = ang_spec_multi_prop(pt, wvl, ...
    delta1, deltan, z, t);

% collimate the beam（去掉输出面上残留的观测面啁啾 e^{+i*pi*r^2/(wvl*R)}，
% 之后相位才有意义——算相干因子必须先准直）
Uvac = Uvac .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));

% 可视化（对照原书图 9.7）
figure;
subplot(1,2,1);
imagesc(xn(1,:)/D2, yn(:,1)/D2, abs(Uvac).^2); axis image; colorbar;
xlabel('2x_n/D_2'); ylabel('2y_n/D_2'); title('真空辐照度 [W/m^2]');
subplot(1,2,2);
imagesc(xn(1,:)/D2, yn(:,1)/D2, angle(Uvac)); axis image; colorbar;
xlabel('2x_n/D_2'); ylabel('2y_n/D_2'); title('真空相位（准直后）[rad]');
