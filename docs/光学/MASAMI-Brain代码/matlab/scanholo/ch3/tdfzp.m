% ===== 书 3.6 节 · Table 3.2 TDFZP.m 的规范化重写：时间依赖菲涅耳波带片的运行条纹 =====
% TDFZP(x,y;z0,t) = A + B*sin[(k0/2*z0)(x^2+y^2) + Omega*t]     (3.6-1)
% 平面波（w0+Omega）与球面波（w0）干涉形成的"会跑的 FZP"：t 增大时同相位圆环半径缩小，
% 条纹向中心运动（running fringes toward the center）——OSH 就是用它光栅扫描 3-D 物体。
% 换瞳（p1=delta, p2=1）得 TDFZP(x,y;z0,-t)：条纹远离中心；sin 全息图不变、cos 全息图反号 (3.6-4)。
%
% 原书参数：B=10.01e6（时间频率，任意）、D=6（空间尺度，任意）、t=linspace(0,1,35)，
% x,y=linspace(-2.5,2.5,256)，FZP = 1+sin(D*(x^2+y^2)+B*t)，逐帧归一化后 movie 播放。
% 本版保留原书参数与 movie 动画，同时输出 3 个静态帧便于笔记截图。
%
% 差异提示：
% 1) 原书三重循环 (ii,jj,kk)；本版用 meshgrid + 广播一次生成（数值相同，速度更快）。
% 2) 原书 image(256*FZP) 行方向向下；本版 imagesc + axis off 展示。
% 3) 改变条纹运动方向：把 sin 里的 +B*t 改成 -B*t（等价于 3.6-4 的 TDFZP(x,y;z0,-t)）。

clear
B = 10.01e6;                 % 时间频率（任意单位）
D = 6;                       % 空间尺度（任意单位）
t = linspace(0, 1, 35);
x = linspace(-2.5, 2.5, 256);
y = linspace(-2.5, 2.5, 256);
[X, Y] = meshgrid(x, y);

% ---- 生成全部帧（原书三重循环的广播版） ----
phase = D*(X.^2 + Y.^2);                 % 空间部分（静态 FZP）
NT = length(t);
FZP = zeros(size(X,1), size(X,2), NT);
for kk = 1:NT
    FZP(:,:,kk) = 1 + sin(phase + B*t(kk));          % TDFZP (3.6-1)
end

% ---- 静态 3 帧（条纹向中心运动：t 增大，同相位圆环半径 r(t)=sqrt((c-B t)/D) 减小） ----
figure(1)
for i = 1:3
    kk = round(1 + (i-1)*(NT-1)/2);
    subplot(1,3,i); imagesc(x, y, 256*FZP(:,:,kk)/max(FZP(:,:,kk),[],'all'));
    colormap(gray(256)); axis image off
    title(sprintf('TDFZP: sin[6r^2 + \\Omega t],  \\Omega t = %.2f', B*t(kk)))
end

% ---- 动画（原书 movie） ----
figure(2)
F = [];
for ll = 1:NT
    subplot(1,1,1); imagesc(x, y, 256*FZP(:,:,ll)); colormap(gray(256)); axis image off
    F = [F, getframe];
end
movie(F, 2)      % 播放 2 遍；条纹向中心收缩（running fringes toward the center）
