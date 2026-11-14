% ===== 书图 7.9 / 7.14:image_fseries f/5、f/10、f/20 系列 =====
% 相干(图 7.9)与非相干(图 7.14)成像的 f/# 系列:理想图像的 fN 固定
% (4.17e5 cyc/m),改出瞳半径只改截止:
%   相干   f0  = wXP/(lam*zXP)   -> f/5: 2e5, f/10: 1e5, f/20: 5e4
%   非相干 2f0 = 1/(lam*(f/#))   -> f/5: 4e5, f/10: 2e5, f/20: 1e5
% 谱显示:相干用 log 幅度(通带外为 0 -> 黑),非相干用五次方根(书口径)。
% 结论(Python 实跑同图样):f/# 越大截止越低、像越糊;同 f/# 下非相干
% 截止是相干的两倍、图样明显更好且无振铃(见笔记 7.2.3/7.3.3 对比)。
% MATLAB/Python 差异:Python 版 ch7_fseries.py;log(0) 在 MATLAB 里是
% -Inf,imagesc 会显示成最小色阶,这里按书做法直接对幅度谱取 log。
% ------------- 脚本骨架(核心行已写全,绘图行可按需展开) -------------
%  image_fseries Coherent/Incoherent image series vs f-number

A=imread('USAF1951_py250','png'); A=flipud(A);
Ig=single(A); Ig=Ig/max(max(Ig));

L=0.3e-3; M=250; du=L/M;
u=-L/2:du:L/2-du;
lambda=0.5*10^-6;
fu=-1/(2*du):1/L:1/(2*du)-(1/L); fv=fu;
[Fu,Fv]=meshgrid(fu,fv);

for fn=[5 10 20]                       %f-number 系列
    f0=1/(2*lambda*fn);                %相干截止(=wXP/(lam*zXP))
    % ---- 相干(图 7.9)----
    H=circ(sqrt(Fu.^2+Fv.^2)/f0);
    ug=sqrt(Ig);
    Hs=fftshift(H); Gg=fft2(fftshift(ug)); Gi=Gg.*Hs;
    ui=ifftshift(ifft2(Gi)); Ii=(abs(ui)).^2;
    % figure: imagesc(u,v,log(abs(fftshift(Gi))+eps)) —— 谱(log)
    % figure: imagesc(u,v,nthroot(Ii,2))           —— 像
    % ---- 非相干(图 7.14)----
    OTF=ifft2(abs(fft2(fftshift(H))).^2); OTF=OTF/OTF(1,1);
    Gg=fft2(fftshift(Ig)); Gi=Gg.*OTF;
    Ii2=ifftshift(ifft2(Gi)); Ii2=real(Ii2); mask=Ii2>=0; Ii2=mask.*Ii2;
    % figure: imagesc(u,v,abs(fftshift(ifft2(fft2(fftshift(Ii2))))).^(1/5))
    % figure: imagesc(u,v,nthroot(Ii2,2))
end
fprintf('fN=4.17e5;f/5:f0=2e5,2f0=4e5;f/10:f0=1e5,2f0=2e5;f/20:f0=5e4,2f0=1e5 (cyc/m)\n');
