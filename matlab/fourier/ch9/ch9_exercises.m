% ch9_exercises.m
% ===== 书 9.4 节习题 9.1-9.6 的 MATLAB 参考实现 =====
% 与 python/fourier/ch9/ch9_exercises.py 数值一致(Python 侧已实跑验证,
% 结果见笔记"习题"小节)。分段复制到命令窗口即可运行。
% 需要:circ.m、rect.m(书附录 C,ch7 目录有 circ.m 副本)、
%       propTF.m(ch5 目录,习题 9.6 用)。
% MATLAB/Python 差异:sinc 归一化定义一致 sin(pi x)/(pi x)(MATLAB 需
%   Signal Processing Toolbox,或用 self-written: sin(pi*x)./(pi*x) 补 x=0);

%% ---- 习题 9.1:矩形线型(式 9.29/9.30) ----
% |gamma(tau)| = |sinc(dnu*tau)| = 0.8270(dnu*tau = 0.3333);
% tau_c = 1/dnu = 0.5ns。把 pc_temp 第 34 行的 S 换成 rect 线型即可。
% 注意:离散谱线权重和 = 1.0196(13 条线落在矩形内,矩形边界在样本之间),
% 仿真对比度比连续理论高 ~1.3%(0.838 vs 0.827)。
lambda0=650e-9; c=3e8; nu0=c/lambda0;
N=51; delnu=2e9; dnu=4*delnu/N;
L1=50e-3; M=250; dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; x1=fftshift(x1); [X1,Y1]=meshgrid(x1,x1);
w=1e-3; dels=5e-3; deld=5e-2; f=0.25; lf=lambda0*f;
tau=deld/c;
fprintf('9.1: |gamma| = |sinc(dnu*tau)| = %.4f; tau_c = 1/dnu = %.2e s\n',...
    abs(sin(pi*delnu*tau)/(pi*delnu*tau)), 1/delnu);
I2=zeros(M);
for n=1:N
   nu=(n-(N+1)/2)*dnu+nu0;
   S=rect((nu-nu0)/delnu)/delnu;          % 式 (9.29):rect 线型
   k=2*pi*nu/c;
   u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)...
       +circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
   u2=1/lf*(fft2(u))*dx1^2;
   I2=I2+S*(abs(u2).^2)*dnu;
end
I2=ifftshift(I2);
x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf;
figure; imagesc(x2,x2,I2); axis square; axis xy; colormap gray;
figure; plot(x2,I2(M/2+1,:)); xlabel('x (m)'); ylabel('Irradiance');

%% ---- 习题 9.2:相干条纹位移 ----
% 相干单频、deld = 0, lam0/4, lam0/2, 3lam0/4:条纹整体平移 deld/lam0 个周期;
% 局限:图样以 lam0 为周期,deld 与 deld+lam0 不可区分(无模糊范围 lam0/2)。
dds=[0 lam0/4 lam0/2 3*lam0/4];
figure;
for jj=1:4
   u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)...
       +circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*2*pi/lam0*dds(jj));
   u2=1/lf*(fft2(u))*dx1^2;
   I=ifftshift(abs(u2).^2);
   subplot(1,4,jj); plot(x2,I(M/2+1,:)); xlim([-1.5e-4 1.5e-4]);
   title(sprintf('\\Deltad = %.2f\\lambda_0', dds(jj)/lam0));
   xlabel('x (m)');
end
fprintf('9.2: 条纹峰随 deld 平移 deld/lam0 个周期;超过 lam0/2 无模糊\n');

%% ---- 习题 9.3:可见度 A1=0.75、A2=0.20(式 9.31) ----
% V = 2*A1*A2/(A1^2+A2^2)*|gamma(deld/c)| = 0.4979*0.6733 = 0.3353。
% (c) 提议解析式:包络用 (1/(lam0*f)^2)*[w^2*J1(z)/(w~r)]^2(无 (9.10) 的
%     双束因子 2),括号换 [A1^2+A2^2+2*A1*A2*|g|cos(...)];A1=A2=1 时退回 (9.10)。
A1=0.75; A2=0.20;
gt=exp(-(pi*delnu*tau/(2*sqrt(log(2))))^2);   % |gamma(deld/c)| = 0.6733
V=2*A1*A2/(A1^2+A2^2)*gt;
fprintf('9.3: V = %.4f\n', V);
I2=zeros(M);
for n=1:N
   nu=(n-(N+1)/2)*dnu+nu0;
   S=1/(sqrt(pi)*b0)*exp(-(nu-nu0)^2/b0^2);  % b0=delnu/(2*sqrt(log(2)))
   k=2*pi*nu/c;
   u=A1*circ(sqrt((X1-dels/2).^2+Y1.^2)/w)...
       +A2*circ(sqrt((X1+dels/2).^2+Y1.^2)/w)*exp(j*k*deld);
   u2=1/lf*(fft2(u))*dx1^2;
   I2=I2+S*(abs(u2).^2)*dnu;
end
I2=ifftshift(I2);
figure; plot(x2,I2(M/2+1,:)); xlim([-2e-4 2e-4]);
xlabel('x (m)'); ylabel('Irradiance'); title(sprintf('V = %.3f', V));

%% ---- 习题 9.4:l_cr = 15/8/4/1 mm ----
% |mu_12| = exp(-(dels/l_cr)^2):0.8948 / 0.6766 / 0.2096 / exp(-25)~0。
% l_cr=15mm 且 mult=2.5 时 sigma_f=37.5mm > L/2=25mm:滤波器通带不足一个
% 频率样本,屏几乎恒定 -> 输出完全相干结果(书描述的失效模式,实测可见度 1.00);
% 调小 21 行倍率(mult 1.5,sigma_f=22.5<25)后可见度降到 ~0.96。
% l_cr=1mm:|mu_12|~0,无条纹,图样退化为双倍单束 Airy 斑;更多屏只更平滑。
lam=650e-9; lf=lam*f;
mults=[2.5 1.5];
for icr=1:2
   Lcr=15e-3; mult=mults(icr);
   sigma_f=mult*Lcr; sigma_r=sqrt(4*pi*sigma_f^4/Lcr^2);
   fx1=-1/(2*dx1):1/L1:1/(2*dx1)-1/L1; fx1=fftshift(fx1);
   [FX1,FY1]=meshgrid(fx1,fx1);
   u1=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)+circ(sqrt((X1+dels/2).^2+Y1.^2)/w);
   F=exp(-pi^2*sigma_f^2*(FX1.^2+FY1.^2));
   I2=zeros(M);
   for n=1:400/2
      fie=(ifft2(F.*(randn(M)+j*randn(M)))*sigma_r/dfx1)*M^2*dfx1^2;
      u2=1/lf*(fft2(u1.*exp(j*real(fie))))*dx1^2; I2=I2+abs(u2).^2;
      u2=1/lf*(fft2(u1.*exp(j*imag(fie))))*dx1^2; I2=I2+abs(u2).^2;
   end
   I2=ifftshift(I2)/400;
   fprintf('9.4: l_cr=15mm mult=%.1f -> 可见度(中心条纹)看剖面图\n', mult);
   figure; plot(x2,I2(M/2+1,:)); xlim([-2e-4 2e-4]);
   title(sprintf('l_{cr}=15mm, mult=%.1f', mult));
end

%% ---- 习题 9.5:单个 5mm x 2mm 矩形孔,l_cr = 20 -> 1 mm ----
% l_cr >> 孔径:相干 sinc x sinc 图样带旁瓣;l_cr 减小:旁瓣先消失、
% 中央瓣变滑;l_cr=1mm << 孔径:几乎非相干,只剩光滑包络。
Lcrs=[20e-3 8e-3 4e-3 1e-3];
figure;
for jj=1:4
   Lcr=Lcrs(jj);
   sigma_f=2.5*Lcr; sigma_r=sqrt(4*pi*sigma_f^4/Lcr^2);
   fx1=-1/(2*dx1):1/L1:1/(2*dx1)-1/L1; fx1=fftshift(fx1);
   [FX1,FY1]=meshgrid(fx1,fx1);
   u1=rect(X1/5e-3).*rect(Y1/2e-3);        % 5mm x 2mm 矩形孔
   F=exp(-pi^2*sigma_f^2*(FX1.^2+FY1.^2));
   I2=zeros(M);
   for n=1:200/2
      fie=(ifft2(F.*(randn(M)+j*randn(M)))*sigma_r/dfx1)*M^2*dfx1^2;
      u2=1/lf*(fft2(u1.*exp(j*real(fie))))*dx1^2; I2=I2+abs(u2).^2;
      u2=1/lf*(fft2(u1.*exp(j*imag(fie))))*dx1^2; I2=I2+abs(u2).^2;
   end
   I2=ifftshift(I2)/200;
   subplot(1,4,jj); imagesc(x2,x2,sqrt(I2)); axis square; axis xy;
   colormap gray; title(sprintf('l_{cr}=%.0fmm', Lcr*1e3));
end

%% ---- 习题 9.6:高斯谢尔模型(GSM)光束 ----
% w0=2cm, lam=1.06um, z=1km, l_cr=1cm, M=256。L1 取 0.55m 满足 TF 判据
% dx >= lam*z/L(2.148mm >= 1.927mm,书 5.4 节),且 sigma_f=2.5cm < L/2。
% 式 (9.33):Delta = 2.721;式 (9.32):I2 = (I0/Delta^2) exp(-2r^2/(w0*Delta)^2),
% I0 = 1(源面峰值辐照度;总功率 pi*w0^2/2 守恒)。
% Python 实测:拟合 Delta_sim = 2.726(理论 2.721);GSM 峰值 0.143 vs 0.135;
% 相干峰值 0.5843 vs 0.5843(逐位一致);总功率 6.2832e-4 m^2 精确守恒。
Mg=256; w0=0.02; lamg=1.06e-6; zg=1000; Lcrg=0.01; L1g=0.55;
k=2*pi/lamg;
xg=-L1g/2:(L1g/Mg):L1g/2-L1g/Mg; [Xg,Yg]=meshgrid(xg,xg);
U1=exp(-(Xg.^2+Yg.^2)/w0^2);                 % 式 (9.34),居中排列
Delta=sqrt(1+(2*zg/(k*w0^2))^2*(1+2*w0^2/Lcrg^2));   % 式 (9.33)
Deltac=sqrt(1+(2*zg/(k*w0^2))^2);
fprintf('9.6: Delta = %.4f(相干 %.4f);峰值 I0/Delta^2 = %.4f / %.4f\n',...
    Delta, Deltac, 1/Delta^2, 1/Deltac^2);
sigma_f=2.5*Lcrg; sigma_r=sqrt(4*pi*sigma_f^4/Lcrg^2); dfx1g=1/L1g;
fxg=-1/(2*(L1g/Mg)):1/L1g:1/(2*(L1g/Mg))-1/L1g; fxg=fftshift(fxg);
[FXg,FYg]=meshgrid(fxg,fxg);
Fg=exp(-pi^2*sigma_f^2*(FXg.^2+FYg.^2));
I2g=zeros(Mg);
for n=1:400/2
   fie=(ifft2(Fg.*(randn(Mg)+j*randn(Mg)))*sigma_r/dfx1g)*Mg^2*dfx1g^2;
   % 注意:propTF 约定源场"居中排列"(内部自己 fftshift),而 fie 是 ifft2
   % 输出(FFT 排列)——相乘前先 ifftshift 回居中,两套排列别混。
   fie=ifftshift(fie);
   u2=propTF(U1.*exp(j*real(fie)),L1g,lamg,zg); I2g=I2g+abs(u2).^2;
   u2=propTF(U1.*exp(j*imag(fie)),L1g,lamg,zg); I2g=I2g+abs(u2).^2;
end
I2g=I2g/400;
u2c=propTF(U1,L1g,lamg,zg);                  % 完全相干对比
figure; imagesc(xg,xg,sqrt(I2g)); axis square; axis xy; colormap gray;
figure; plot(xg,I2g(Mg/2+1,:),xg,1/Delta^2*exp(-2*xg.^2/(w0*Delta)^2),'--',...
    xg,abs(u2c(Mg/2+1,:)).^2,':');
xlabel('x (m)'); legend('GSM sim','GSM analytic (9.32)','coherent');
fprintf('9.6: GSM 总功率 %.4e m^2(源 pi*w0^2/2 = %.4e,守恒)\n',...
    sum(sum(I2g))*(L1g/Mg)^2, pi*w0^2/2);
