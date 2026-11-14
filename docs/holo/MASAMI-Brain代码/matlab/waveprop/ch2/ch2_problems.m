%% ===== 书 Schmidt PM199 第 2 章 =====
%  原书 2.7 节习题 1-5：五个信号在 a=1 与 a=10 两档下的 DFT vs 解析 FT
%  与 python/waveprop/ch2/ch2_problems.py 一致（数值以 Python 实跑为准）。
%  采样策略（原书 2.5 节案例研究的模板）：
%    sinc：奈奎斯特 delta=0.85/a；高斯类：p=0.01 电平带宽；
%    tri / exp(-a|x|)：谱 1/f^2 衰减，取 |G|=0.01 峰值为 f_max。
%  需要 R2016b+。
clear; clc; close all;

a_list = [1, 10];
p = 0.01; b3 = 2;
names = {'1: sinc(ax)', '2: exp(-pi a^2 x^2)', '3: quad-phase, b=2', ...
         '4: tri(ax)', '5: exp(-a|x|)'};
fprintf('===== 习题 1-5：a=1 / a=10 的采样参数与 DFT 精度 =====\n');
fprintf('%-18s %-5s %-10s %-8s %-12s %-12s\n', ...
    '信号', 'a', 'delta', 'N', '带内误差', '带外误差');
rows = zeros(10, 8);   % [a delta N err_in err_out ...] 按题号-a 排列
k = 0;
for q = 1:5
    for a = a_list
        switch q
            case 1,  delta = 0.85/a;                       % 奈奎斯特 + 15%
            case 4,  f_max = a/(pi*sqrt(p));   delta = 1/(2*f_max);
            case 5,  f_max = a*sqrt(1/p - 1)/(2*pi); delta = 1/(2*f_max);
            otherwise
                if q == 2, coeff = a^2; else, coeff = a^2 + b3^4/a^2; end
                f_max = sqrt(coeff/pi*log(1/p)); delta = 1/(2*f_max);
        end
        switch q
            case 4, half = 1.2/a;
            case 1, half = 12/a;
            case 5, half = log(100)/a;
            otherwise, half = sqrt(log(100)/pi)/a;
        end
        N = 2^ceil(log2(2*half/delta));                % 取 2 的幂
        [x, f] = coords(N*delta, N);
        g = signal_eval(q, x, a, b3);
        G = ft(g, delta);
        G_a = analytic_eval(q, f, a, b3);
        dev = abs(abs(G) - abs(G_a)) / max(abs(G_a));
        if q == 1
            inb = abs(f) < a/2;
            e_in = max(dev(inb)); e_out = max(dev(~inb));
        else
            e_in = max(dev); e_out = NaN;
        end
        k = k + 1;
        rows(k,:) = [q, a, delta, N, e_in, e_out, 0, 0];
        fprintf('%-18s %-5.0f %-10.4g %-8d %-12.2e %-12.4g\n', ...
            names{q}, a, delta, N, e_in, e_out);
    end
end

%% 绘图：5 行 x 2 列
figure('Position',[30 30 1150 1500]);
k = 0;
for q = 1:5
    for c = 1:2
        k = k + 1;
        a = a_list(c);
        delta = rows(k,3); N = rows(k,4);
        [x, f] = coords(N*delta, N);
        G = abs(ft(signal_eval(q, x, a, b3), delta));
        G_a = abs(analytic_eval(q, f, a, b3));
        subplot(5,2,k);
        m = abs(f) <= 0.8/delta;
        plot(f(m), G_a(m), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); hold on;
        plot(f(m), G(m), 'b.', 'MarkerSize', 4);
        grid on;
        title(sprintf('%s  a=%g  (delta=%.3g, N=%d)', names{q}, a, delta, N), ...
            'FontSize', 10);
        if q == 5, xlabel('f_x [m^{-1}]'); end
        if q == 1 && c == 1
            legend({'解析', 'DFT'}, 'FontSize', 9, 'Location', 'south');
        end
    end
end

%% ===== 本地函数 =====
function G = ft(g, delta)
    G = fftshift(fft(fftshift(g))) * delta;            % 原书 Listing 2.1
end

function [x, f] = coords(L, N)
    delta = L / N;
    n = (-N/2 : N/2-1);
    x = n * delta;
    f = n / (N * delta);
end

function s = mysinc(x)
    s = ones(size(x));
    nz = x ~= 0;
    s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function y = tri_fun(x)
    % 三角函数 Lambda(x)：|x|<=1 取 1-|x|，其余 0（原书附录 A）
    y = double(abs(x) <= 1) .* (1 - abs(x));
end

function g = signal_eval(q, x, a, b)
    switch q
        case 1,  g = mysinc(a*x);
        case 2,  g = exp(-pi*(a*x).^2);
        case 3,  g = exp(-pi*(a*x).^2) .* exp(1i*pi*(b*x).^2);
        case 4,  g = tri_fun(a*x);
        case 5,  g = exp(-a*abs(x));
    end
end

function G = analytic_eval(q, f, a, b)
    switch q
        case 1,  G = (abs(f) < a/2)/a + (abs(f) == a/2)*(0.5/a);
        case 2,  G = 1/a * exp(-pi*(f/a).^2);
        case 3,  c = a^2 - 1i*b^2;
                 G = 1/sqrt(c) * exp(-pi*f.^2/c);
        case 4,  G = mysinc(f/a).^2 / a;                   % tri 的 FT：sinc^2/a
        case 5,  G = 2*a ./ (a^2 + (2*pi*f).^2);
    end
end
