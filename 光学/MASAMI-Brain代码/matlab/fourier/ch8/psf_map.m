% psf_map generate psf map
% ===== 书 8.7.1 节:psf_map 像面 PSF 分布图 =====
% 书 8.7.1 节列出的脚本(原书第 1-53 行)。f/20 出瞳,M=250、L=1 mm,
% 2f0 = 1e5 1/m <= fN = 1.25e5 1/m,式 (8.6) 满足(书评注:分号故意不写,
% 便于改参数时核对两个频率)。像差夸大:0.5λ 球差 + 1λ 彗差 + 1.5λ 像散。
% 思路:对 7x7 个归一化像点 (u0,v0) 各生成 PSF(式 8.7),circshift 平移到
% 像面对应位置后叠加显示。注意书 8.7.1 节评注 (f):这只是"摆放",不是卷积
% 或叠加运算(叠加见 8.7.2 的 image_super)。
% 省时技巧(评注 b):第 26 行先把频率坐标 fu=fftshift(fu) 预移位,循环内
% H 天然处于 FFT 序,第 41 行 ifft2 前不再需要 fftshift。
% circshift 评注 (e):第一个变量是行(= y/v 方向),round 取整;移出边界的
% 元素从对侧回来(周期延拓概念)。
% MATLAB/Python 差异:Python 版 python/fourier/ch8/ch8_psf_map.py 用
%   np.roll(h2, (r, c), axis=(0,1)),轴含义与 circshift 的 [行,列] 相同;
%   numpy 的 round 是银行家舍入,本题 ±87.5 恰好也舍入到 ±88,结果一致。
% ------------- 以下为书原代码(第 1-53 行) -------------
% psf_map generate psf map

M=250;
L=1e-3;                %image plane side length
du=L/M;                %sample interval
u=-L/2:du:L/2-du; v=u; %coordinates

lambda=0.5*10^-6;      %wavelength
k=2*pi/lambda;         %wavenumber
wxp=2.5e-3;            %exit pupil radius
zxp=100e-3;            %exit pupil distance
fnum=zxp/(2*wxp)       %exit pupil f-number

twof0=1/(lambda*fnum)  %inc cutoff freq
fN=1/(2*du)            %Nyquist frequency

% aberration coefficients
wd=0*lambda;
w040=0.5*lambda;
w131=1*lambda;
w222=1.5*lambda;
w220=0*lambda;
w311=0*lambda;

fu=-1/(2*du):1/L:1/(2*du)-(1/L); %image freq coords
fu=fftshift(fu); %shift cords, avoid shift H in loop
[Fu,Fv]=meshgrid(fu,fu);

I=zeros(M);
% loop through image plane positions
for u0=[-.7:.7/3:.7]
    for v0=[-.7:.7/3:.7]
        % wavefront
        W=seidel_5(u0,v0,-2*lambda*fnum*Fu...
            ,-2*lambda*fnum*Fv,...
            wd,w040,w131,w222,w220,w311);
        % coherent transfer function
        H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum)...
            .*exp(-j*k*W);
        % PSF
        h2=abs(ifftshift(ifft2(H))).^2;
        % shift PSF to image plane position
        h2=circshift(h2,[round(v0*M/2)...
            ,round(u0*M/2)]);
        % add into combined frame
        I=I+h2;
    end
end

figure(1)
imagesc(u,v,nthroot(I,4));
xlabel('u (m)'); ylabel('v (m)');
colormap('gray'); axis square; axis xy
