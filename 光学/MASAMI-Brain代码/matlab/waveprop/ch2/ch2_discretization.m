%% ===== 书 Schmidt PM199 第 2 章 =====
%  离散化的三种失真与三条补救路线（对应原书 2.3、2.4 节与图 2.6-2.9）
%  信号 g(x)=exp(-a|x|)（a=10），G(f)=2a/(a^2+(2*pi*f)^2)。
%  链条：采样(混叠) -> 截断(波纹) -> 频域采样(空域周期延拓)。
%  中间谱的精确算法：S(f) = sum_n g(n*delta) exp(-i2*pi*f*n*delta)，
%  DFT 输出 = delta*S(f_m)（脚本内自检）。
%  与 python/waveprop/ch2/ch2_discretization.py 一致；需要 R2016b+。
clear; clc; close all;

a = 10.0;
strat = {  % 名称, delta, L, N
    '基线（原书图 2.6）',              0.0375,  0.6, 16
    'A 加大网格：\delta x2（图 2.7）',  0.075,   1.2, 16
    'B 加密采样：N x2（图 2.8）',       0.01875, 0.6, 32
    'C 双管齐下（图 2.9）',             0.01875, 1.2, 64
};

%% 链条自检 + 误差量化
f_fine = linspace(-20, 20, 8001);
G0 = 2*a / a^2;                          % 解析峰值 2/a = 0.2
n_ali = -16:15; x_ali = n_ali/8;
fprintf('===== 混叠折返恒等式（f1=6, delta=1/8 m -> f2=2 m^-1）=====\n');
fprintf('    max|cos(2*pi*6*x) 与 cos(2*pi*2*x) 的样本差| = %.1e\n', ...
    max(abs(cos(2*pi*6*x_ali) - cos(2*pi*2*x_ali))));

fprintf('===== 三种失真的量化（频带 |f| <= 1/(2 delta)，误差/解析峰值）=====\n');
fprintf('%-28s %-12s %-12s %-12s\n', '策略', '混叠误差', '波纹误差', 'DFT总误差');
results = cell(size(strat,1), 1);
for s = 1:size(strat,1)
    name = strat{s,1}; delta = strat{s,2}; L = strat{s,3}; N = strat{s,4};
    n = (-N/2 : N/2-1);
    x_n = n * delta;
    g_samp = exp(-a*abs(x_n));                           % 居中排列
    % 中间谱 S(f) = sum_n g(n*delta) exp(-i2 pi f n delta)（矩阵一行一个频率）
    S = real(exp(-2i*pi*x_n.' * reshape(f_fine,1,[])) * g_samp.');  % (1 x Mf)
    G_ali = alias_only(f_fine, a, delta);                % (1/delta)*sum_k G
    Gt = G_true(f_fine, a);
    fmax = 1/(2*delta);
    band = abs(f_fine) <= fmax;
    e_alias = max(abs(delta*G_ali(band) - Gt(band))) / G0;
    e_rip   = max(abs(delta*(S(band) - G_ali(band)))) / G0;
    fm = (-N/2 : N/2-1) / (N*delta);
    G_dft = ft(g_samp, delta);
    e_dft = max(abs(G_dft - G_true(fm,a))) / G0;
    fprintf('%-28s %-12.2e %-12.2e %-12.2e\n', name, e_alias, e_rip, e_dft);
    results{s} = struct('name',name,'delta',delta,'L',L,'N',N, ...
        'f',f_fine,'Gt',Gt,'Ga',G_ali,'S',S, ...
        'fm',fm,'Gdft',abs(G_dft),'e',e_dft);
end

%% 绘图：4 行 x 2 列（左列 1/delta 放大标度，同原书图 (d)(f)）
figure('Position',[40 40 1150 1250]);
for s = 1:4
    r = results{s};
    subplot(4,2,2*s-1);
    plot(r.f, r.Gt/r.delta, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.4); hold on;
    plot(r.f, r.Ga, ':', 'Color', [0.2 0.4 0.85], 'LineWidth', 1.0);
    plot(r.f, r.S, 'k-', 'LineWidth', 1.1);
    xlim([-20 20]); ylim([-0.05 1.12*max(r.Ga)]);
    title(sprintf('%s   \\delta=%.4g, L=%.2g, N=%d', r.name, r.delta, r.L, r.N), ...
        'FontSize', 10);
    ylabel('幅值（1/\delta 放大标度）'); grid on;
    if s == 1
        legend({'解析 G(f_x)/\delta', 'G\otimescomb(\delta f_x)（混叠）', ...
                '再\otimesL sinc(L f_x)（+波纹）'}, 'FontSize', 8, 'Location', 'northeast');
    end
    if s == 4, xlabel('f_x [m^{-1}]'); end

    subplot(4,2,2*s);
    stem(r.fm, r.Gdft, 'b.', 'MarkerSize', 4); hold on;
    plot(r.f, r.Gt, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.4);
    xlim([-20 20]); ylim([-0.05 0.32]); grid on;
    title(sprintf('DFT 输出 vs 解析（总误差 %.1e）', r.e), 'FontSize', 10);
    if s == 1, legend({'DFT', '解析'}, 'FontSize', 8, 'Location', 'northeast'); end
    if s == 4, xlabel('f_x [m^{-1}]'); end
end

%% ===== 本地函数 =====
function G = G_true(f, a)
    % 解析 FT：G = 2a/(a^2+(2 pi f)^2)（原书式 2.14）
    G = 2*a ./ (a^2 + (2*pi*f).^2);
end

function Ga = alias_only(f, a, delta)
    % 只做第 1 步的频谱：(1/delta) * sum_k G(f + k/delta)
    Ga = zeros(size(f));
    for k = -400:400
        Ga = Ga + G_true(f + k/delta, a);
    end
    Ga = Ga / delta;
end

function G = ft(g, delta)
    % 一维 DFT（原书 Listing 2.1）
    G = fftshift(fft(fftshift(g))) * delta;
end
