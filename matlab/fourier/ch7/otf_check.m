% ===== otf_check:OTF 的 FFT 自相关 vs 解析式(7.29)/(7.31)+ 免缠绕条件 =====
% 笔记 7.3 节"拓展:式 (7.17) 的第二张面孔"的数值验证(Python 实跑结果):
%   OTF 用 FFT 算 = 循环自相关 = 真自相关的周期延拓(副本间距 1/du);
%   免缠绕条件 2*f0 <= 1/du - 2*f0,即 2*f0 <= fN —— 恰为判据 (7.17)。
% 实测(Python):f/10(2f0=2e5 <= fN=4.17e5):圆瞳 OTF 与式 (7.31) 最大差
%   3.3e-3、方瞳与式 (7.29) 最大差 1.6e-2;改 f/4(2f0=5e5 > fN)最大差
%   7.9e-2 —— 高频区被环绕污染。OTF 支撑半径实测恰为 2f0。
% MATLAB/Python 差异:Python 版 ch7_incoh_image.py;jinc 类解析式里
%   arccos 的宗量要 clip 到 [0,1]。
% ------------- 脚本骨架 -------------
%  otf_check  OTF analytic vs FFT-autocorrelation

L=0.3e-3; M=250; du=L/M;
lambda=0.5*10^-6;
fu=-1/(2*du):1/L:1/(2*du)-(1/L);
[Fu,Fv]=meshgrid(fu,fv);

for fn=[10 4]                          %f/10 合格;f/4 违反式 (7.17)
    f0=1/(2*lambda*fn);                %相干截止
    H=circ(sqrt(Fu.^2+Fv.^2)/f0);
    OTF=ifft2(abs(fft2(fftshift(H))).^2); OTF=OTF/OTF(1,1);
    OTFc=fftshift(real(OTF));          %转居中排列

    % 式 (7.31) 圆瞳解析 OTF:rho = sqrt(Fu^2+Fv^2),r = rho/(2*f0)
    rho=sqrt(Fu.^2+Fv.^2);  r=min(rho/(2*f0),1);
    OTFa=(2/pi)*(acos(r)-r.*sqrt(1-r.^2)); OTFa(rho>2*f0)=0;

    err=max(max(abs(OTFc-OTFa)));
    fprintf('f/%g: 2f0=%.3g, fN=%.3g, max|FFT-解析| = %.2e\n', ...
            fn, 2*f0, 1/(2*du), err);
    % figure: plot(fu,OTFc(M/2+1,:),fu,OTFa(M/2+1,:),'--')  沿 fV=0 切片
end

% 方瞳对照(式 7.29):H=rect(fU/2f0)rect(fV/2f0) -> OTF=Lambda(fU/2f0)Lambda(fV/2f0)
f0=1/(2*lambda*10);
rect_a=@(x) double(abs(x)<=1/2);       %rect(基础 MATLAB 没有,附录 C.3)
Hs=rect_a(Fu/(2*f0)).*rect_a(Fv/(2*f0));
OTFs=ifft2(abs(fft2(fftshift(Hs))).^2); OTFs=OTFs/OTFs(1,1);
Lam=@(t)(abs(t)<=1).*(1-abs(t));       %Lambda 匿名函数
OTFsa=Lam(Fu/(2*f0)).*Lam(Fv/(2*f0));
err_s=max(max(abs(fftshift(real(OTFs))-OTFsa)));
fprintf('方瞳:f/10, max|FFT-解析(7.29)| = %.2e(Python 实测 1.6e-2,rect 边缘 Gibbs)\n', ...
        err_s);
