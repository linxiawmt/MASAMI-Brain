% ch8_exercises.m
% ===== 书 8.8 节习题 8.2/8.3/8.5/8.8 的 MATLAB 参考实现 =====
% 与 python/fourier/ch8/ch8_exercises.py 数值一致(Python 侧已实跑验证,
% 结果见笔记"习题"小节)。分段复制到命令窗口即可运行。
% 需要:seidel_5.m(本目录)、circ.m/rect.m(书附录 C)。

%% ---- 习题 8.2:W040=1λ 时用离焦补偿球差 ----
% 理论:var(c) = c^2/12 + c/6 + 4/45(c = Wd/W040)在 c = -1 最小,
% 即平衡形式 W040(r^4 - r^2),σ_min = W040/sqrt(180) ≈ 0.0745*W040。
M=512; L=1e-3; lambda=0.55e-6; k=2*pi/lambda;
wxp=0.01; zxp=0.1; lz=lambda*zxp; f0=wxp/(lambda*zxp);
x=-1:2/(M-1):1; [X,Y]=meshgrid(x,x);
inside=(sqrt(X.^2+Y.^2)<=1);
w040=1*lambda;
Wd_list=-1.5:0.005:0;
rms=zeros(size(Wd_list));
for i=1:numel(Wd_list)
    W=seidel_5(0,0,X,Y,Wd_list(i)*lambda,w040,0,0,0,0);
    rms(i)=sqrt(mean(W(inside).^2)-mean(W(inside))^2);
end
[~,imin]=min(rms);
fprintf('8.2: 最优 Wd = %+.3f lam, min rms = %.4f lam (理论 -1 与 1/sqrt(180))\n',...
    Wd_list(imin)/lambda, rms(imin)/lambda);

%% ---- 习题 8.3:through-focus PSF(δd = 8(f/#)^2 * Wd)----
% f/5 透镜、W040=4.963λ:δd=-0.5mm(Wd=-4.55λ)最接近平衡点
% δd = 8(f/#)^2*(-W040) = -0.55mm,PSF 峰值最高(约为 δd=0 的 5 倍)。
M=1024; L=1e-3; lambda=0.55e-6; k=2*pi/lambda;
wxp=0.01; zxp=0.1; lz=lambda*zxp; f0=wxp/(lambda*zxp); fnum=5;
du=L/M; u=-L/2:du:L/2-du;
fu=-1/(2*du):1/L:1/(2*du)-1/L; [Fu,Fv]=meshgrid(fu,fu);
Xh=-lz*Fu/wxp; Yh=-lz*Fv/wxp;
dds=[-1 -0.5 0 0.5 1];                 % mm
for i=1:numel(dds)
    Wd=dds(i)*1e-3/(8*fnum^2);
    W=seidel_5(0,0,Xh,Yh,Wd,4.963*lambda,0,0,0,0);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    fprintf('8.3: dd=%+.1f mm -> peak = %.4g\n', dds(i), max(h2(:)));
end

%% ---- 习题 8.5:最大像高点的临界采样离焦 ----
% 式 (8.12) 界 L/(4f#) = 50um = 90.91λ;式 (8.16) 现有斜率(û0=1)= 61.04λ
% -> 临界 Wd = (90.91-61.04)/2 = +14.93λ;超过后波前相位欠采样,PSF 出伪像。
cf=[0 4.963 2.637 9.025 7.536 0.157]*lambda;   % [wd w040 w131 w222 w220 w311]
bound=L/(4*fnum);
slope0=4*cf(2)+3*cf(3)+2*cf(4)+2*cf(5)+cf(6);  % û0=1、x̂=1,ŷ=0
wd_crit=(bound-slope0)/2;
fprintf('8.5: 临界 Wd = %+.2f lam\n', wd_crit/lambda);
for fac=[0.5 1 1.5]
    cft=cf; cft(1)=fac*wd_crit;
    W=seidel_5(0,1,Xh,Yh,cft(1),cft(2),cft(3),cft(4),cft(5),cft(6));
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    fprintf('      Wd x%.1f -> peak = %.4g\n', fac, max(h2(:)));
end

%% ---- 习题 8.8:Strehl 比值 ----
% (a) 无像差峰值 Iunab;(b) W040=0..1λ 的峰值;(c) S 曲线;
% (d) Marechal 近似 S≈exp(-(kσ)^2),σ^2 = (1/12+1/180)*W040^2
%     (书印本漏了 W040 上的平方,量纲分析见笔记勘误 4);
% (e) f/5 例镜轴上:S ≈ 0.021(Python 实测;Marechal 在 5λ 处失效)。
H0=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp);
h2unab=abs(ifftshift(ifft2(fftshift(H0)))).^2;
Iunab=max(h2unab(:));
wvs=[0 0.1 0.25 0.5 0.75 1.0];
S=zeros(size(wvs));
for i=1:numel(wvs)
    W=seidel_5(0,0,Xh,Yh,wvs(i)*lambda,0,0,0,0,0);
    H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
    h2=abs(ifftshift(ifft2(fftshift(H)))).^2;
    S(i)=max(h2(:))/Iunab;
end
S_marechal=exp(-(k^2)*(1/12+1/180)*(wvs*lambda).^2);
fprintf('8.8:  W040(lam)   S_num    S_marechal\n');
for i=1:numel(wvs)
    fprintf('      %5.2f     %.4f    %.4f\n', wvs(i), S(i), S_marechal(i));
end
figure; semilogy(wvs,S,'o',wvs,S_marechal,'-');
xlabel('W_{040} (waves)'); ylabel('Strehl ratio');
legend('numerical','Marechal'); grid on
