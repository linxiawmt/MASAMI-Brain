function fresnel_lsi
% ===== 拓展 1.2.3 · 自由空间当作 LSI 系统：菲涅耳衍射的频域实现 =====
% 改编自你已有的 Fresnel Diffraction.m（规范化重写，物理参数与结果一致）：
% 矩形孔 f(x,y) 放在 z=0，求 z=5 cm 处的衍射光强。
% 本章视角（1.2 节）：自由空间是一个 LSI 系统，
%   空间域：g = f * h，h 为菲涅耳脉冲响应（二次相位 + 高斯窗）；
%   频域：  G = F*H（卷积定理）——本脚本按原书流程用 FFT 实现。
% MATLAB / Python 差异提示：本脚本与 fresnel_lsi.py 逐行对应；
%   meshgrid 后 h(r,c) 里用 X(r),Y(c)（与原书一致），Python 侧用广播实现。

L=1; % 显示区域 L x L（单位 cm）
N=256; % 采样点数
dx=L/(N-1); % 采样步长 1/255 cm
M=111; % 正方形孔的像素宽度（111*dx 约为 0.43 cm）
lambda0=0.00006328; % He-Ne 红光波长（cm），即 632.8 nm
z=5; % 传播距离（cm）

% ---- 输入：中心 M x M 正方形孔（与书中 fft2Drect.m 同样的做法） ----
R=zeros(256);
r=ones(M);
n=(M-1)/2;
R(128-n:128+n,128-n:128+n)=r;

% ---- 坐标轴标定 ----
for k=1:256
    X(k)=1/255*(k-1)-L/2;   % 空间坐标 [-0.5, 0.5] cm
    Y(k)=1/255*(k-1)-L/2;
end

% ---- 菲涅耳脉冲响应 h(x,y;z)：二次相位 + 高斯窗 ----
sigma=pi/(lambda0*z);       % sigma = pi/(lambda*z)，二次相位系数
for rr=1:256
    for cc=1:256
        % j*(sigma/pi) = j/(lambda*z) 振幅因子；
        % exp(-4*200*(x^2+y^2)) 是高斯低通窗（原书写法），把脉冲响应截断在
        % 显示窗口内，抑制离散卷积的卷绕误差。
        h(rr,cc)=j*(sigma/pi)*exp(-4*200*(X(rr).^2+Y(cc).^2))* ...
                 exp(-j*sigma*(X(rr).^2+Y(cc).^2));
    end
end

% ---- 频域实现卷积（原书流程：FR.*H 后 ifft2 再 fftshift） ----
FR=(1/256)^2*fft2(R);
H=(1/256)^2*fft2(h);
HR=FR.*H;
g=ifft2(HR);
g=(256^2)*g;
g=fftshift(g);
I=abs(g).^2;                % 光强

% ---- 可视化（对应原书 figure(1)(2)(3) 与 CCD 视图） ----
figure(1)
image(X,Y,255*R);
colormap(gray(256)); axis square; axis image
xlabel('cm'); ylabel('cm')
title('Input square aperture rect(x/a, y/a)')

figure(2)
plot(X+dx/2,R(:,127))
grid; axis([-0.5 0.5 -0.1 1.2])
xlabel('cm')
title('rect function: cross-section plot')

figure(3)
plot(X+dx/2,abs(g(:,127)))
grid
xlabel('cm')
title('|g(x,y)| cross-section at z = 5 cm')

figure(4)
imagesc(X,Y,I); axis image
colormap(gray(256)); colorbar
xlabel('x (cm)'); ylabel('y (cm)')
title('2-D diffraction pattern (CCD view) at z = 5 cm')

end
