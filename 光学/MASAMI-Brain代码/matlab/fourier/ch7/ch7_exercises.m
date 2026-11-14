% ===== 书 7.4 节习题 7.1-7.10:MATLAB 骨架 =====
% 数值结果与图见笔记正文与 python/fourier/ch7/ch7_exercises.py 实跑。
% 依赖:circ.m、udelta.m(本目录)、make_usaf_chart.m 生成的 PNG。

%% ---- 习题 7.1:薄透镜成像(r=10 mm,f=50 mm,物距 200 mm)----
% (a) 透镜定律(式 7.1):1/f = 1/z1 + 1/z2
%     1/z2 = 1/50 - 1/200 = 3/200 -> z2 = 66.67 mm
%     横向放大率(式 7.2):Mt = -z2/z1 = -1/3(倒像,缩小)
% (b) 薄透镜:瞳=透镜 -> DXP = 20 mm、zXP = z2 = 66.67 mm;
%     f/# = zXP/DXP = 66.67/20 = 3.33
% (c) z1 >> f:1/z2 = 1/f - 1/z1 -> z2 -> f = 50 mm(像聚到焦面)

%% ---- 习题 7.2:DXP=1/2 in、zXP=4 in、lam=0.587 um ----
% (a) f/# = 4/0.5 = 8;f0 = wXP/(lam*zXP) = 1.06e5 cyc/m;
%     2f0 = 1/(lam*(f/#)) = 2.13e5 cyc/m
% (b) du <= lam*(f/#)/2 = 2.348 um(式 7.18)
% (c) L = 1 mm -> M >= L/du = 426(实用取 512)

%% ---- 习题 7.3:矩形出瞳(wx=1 mm、wy=0.5 mm、zXP=50 mm、L=1 mm、M=250)----
% (a) 截止:f0x = wx/(lam*zXP) = 4e4、f0y = wy/(lam*zXP) = 2e4;
%     非相干 2f0x = 8e4、2f0y = 4e4 cyc/m
% (b) 相干:du=4 um;H=rect(Fu/2f0x).*rect(Fv/2f0y)(基础 MATLAB 无 rect,
%     可用逻辑掩码:mask=(abs(Fu)<=2*f0x)&(abs(Fv)<=2*f0y));
%     Gg=fft2(fftshift(ug)); Gi=Gg.*fftshift(H); ui=ifftshift(ifft2(Gi));
%     surf(fu,fv,H) + imagesc(u,v,sqrt(abs(ui).^2))
% (c) 非相干:OTF=ifft2(abs(fft2(fftshift(H))).^2); OTF=OTF/OTF(1,1);
%     Ii=ifftshift(ifft2(fft2(fftshift(Ig)).*OTF)); 谱支撑实测恰为 8e4/4e4
% (d) 分辨率各向异性:竖条(u 向频率,截止 4e4)清楚、横向细节(v 向,
%     截止 2e4)被抹掉——本板组-1/单元2 = 3.8e4 cyc/m:u 向可分辨、
%     v 向不可分辨,一眼可见(Python 图 out_fig7_ex73.png)
% 坑:书 coh_image 第 34 行 xlabel('fv')/ylabel('fu') 与 surf 映射相反,
%     圆瞳无影响、矩形瞳必须改正!

%% ---- 习题 7.4:环形出瞳(wo=1 mm、wi=0.5 mm、zXP=50 mm)----
% (a) 相干截止由外半径决定:f0 = wo/(lam*zXP) = 4e4(2f0 = 8e4);
%     内半径挖掉 |f| < wi/(lam*zXP) = 2e4 的 DC 与低频
% (b) 相干:H = (rho<4e4) - (rho<2e4)(环带);
%     图像边缘增强(伪相衬):大条内部平均场 0.528 -> 0.101(实测)
% (c) 非相干:OTF = 环的自相关归一 —— 低频受损、中高频相对占优
%     (实测与全孔径之比:|f|~3e4 处 0.48、|f|~6e4 处 1.29)
% (d) 丢 DC/低频的效果:平坦区域的场被抽掉,只剩边缘 -> "只有轮廓没有
%     内部";大尺度亮度信息丢失,人眼看是"描边"

%% ---- 习题 7.5:粗糙物体(见 coh_image_rough.m)----
% (a) 相位尺度 2*pi -> pi -> 0.4*pi:斑纹对比度 C=std/mean 断崖下跌
%     (均匀物体实测 0.998/0.20/0.03),测试板重新可辨。
%     解释:各点波前相位差变小,干涉起伏(斑纹)变弱。
% (b) 孔径变小(f/5 -> f/20):PSF 变宽 -> 斑纹颗粒变大,
%     颗粒线度约 1.22*lam*(f/#)(实测 3.6/6.0/9.6 um,网格 1.2 um)。

%% ---- 习题 7.6:衍射限 PSF(Airy)----
% (a) h(u,v) = F^-1{circ(rho/f0)}(式 7.13 的逆变换):
%     h = pi*f0^2 * 2*J1(2*pi*f0*u)/(2*pi*f0*u)  (jinc 型,u=sqrt(u^2+v^2))
%     MATLAB:用 base 思路 h=fftshift(fft2(ifftshift(H)))*(df)^2,H 实偶
%     时 FT=IFT(与 Python 的 base.ft2_centered 同一定标)
% (b) PSF = |h|^2 = (pi*f0^2)^2*[2J1(v)/v]^2,v = 2*pi*f0*u
%     = pi*D*u/(lam*zXP) —— 与式 (4.35) 的 Airy 图样同形
%     (f0 = wXP/(lam*zXP) -> v = 2*pi*wXP*u/(lam*zXP))
% (c) J1 第一零点 v=3.8317 -> u0 = 0.61/f0 = 1.22*lam*(f/#);
%     中央瓣全宽 = 2*u0 = 2.44*lam*(f/#)
% 实测(Python,f/10):h(0) 差 0.65%(DFT 近似),第一零点 6.14 um
% (理论 6.10),PSF vs 解析最大差 2.2e-3。

%% ---- 习题 7.7:Rayleigh 判据(f/# = 20)----
% delta = 1.22*lam*(f/#) = 12.2 um;du = 1.2 um -> S = delta/du = 10.17 -> 10
% Ig=zeros(M); Ig(M/2+1,M/2+1-S/2)=1; Ig(M/2+1,M/2+1+S/2)=1;   %(c)
% wxp = zXP/(2*20) = 3.125 mm -> f0 = 5e4;OTF 同 incoh_image;
% 两点像剖面:中央凹陷/峰 = 0.760(实测)vs Airy 理论 0.735 -> 刚可分辨。

%% ---- 习题 7.9:相衬成像(图 7.15)----
% (a) theta=sqrt(Ig)*pi/100; ug=exp(j*theta);(条形处相位 pi/100,背景 0)
% (b) 直接辐照度 |ug|^2 = 1 —— 常数,肉眼不可见(这就是相衬要解决的问题)
% (c) G=fft2(ug); G(1,1)=G(1,1)*exp(j*pi/2); ui=ifft2(G); Ii=abs(ui).^2;
%     (MATLAB 的 DC 点是 (1,1);numpy 是 [0,0])
%     实测:条纹上 1.059 vs 条纹外 0.996,差 +0.064 ≈ 2*theta = 0.063
%     —— 正相衬 Ii ≈ 1+2*theta
% (d) beta=3*pi/2:条纹上 0.942 vs 条纹外 1.004,差 -0.062 —— 负相衬
% (e) beta=pi:DC 反相 -> 线性衬度消失(只剩 theta^2 二阶小量,max|Ii-1|
%     = 1.6e-3);beta=0:Ii 恒 1。物理:DC(平均场)与相位结构差 pi/2,
%     相移后同相/反相干涉,把相位差转成辐照度差。
% 注意:像面上的平均场修正 —— 精确说是 Ii ≈ |1 + j*theta|^2 型,
%     绝对值含 <ug> 因子,但"条纹-背景差 = 2*theta"与理论一致(实测)。

%% ---- 习题 7.10:傅里叶全息(图 7.16/7.17)----
% (a) M=500、L=10 mm(dx=0.02 mm);物体(F 形,尺寸见书图 7.17,mm):
%     u1 = rect((X+2.0*mm)/0.4mm)*rect(Y/2mm) ...                    %竖梁
%        + rect((X+1.3*mm)/1mm)*rect((Y-0.8*mm)/0.4mm)               %上臂
%        + rect((X+1.3*mm)/1mm)*rect(Y/0.4mm);                       %中臂
%     点源(参考):1000*udelta((X-2*mm)/dx).*udelta(Y/dx)
% (b) 点源乘 1000(亮参考)
% (c) 记录:G=fft2(u1); T=abs(G).^2;(坐标与幅度定标按书忽略)
% (d) 回放:P=ifft2(T);(书用逆 FFT 保持方向;实际回放透镜做正变换,
%     像转 180 度,用 ifft2 抵消)
% (e) P(1,1)=0;  %轴上亮点 = R^2*delta(回放中仍是平面波的部分),显示干扰
% (f) imagesc 显示 fftshift(abs(P))(居中排列)
% 结果(实测):真像在 x = -3.5 mm(物中心 -1.5 - 参考偏移 2;与物体
%     归一化互相关 +1.000),共轭像(180 度翻转)在 +3.5 mm(与翻转物体
%     互相关 +0.929),中间是物体自相关宽包络 —— 傅里叶自相关定理:
%     ifft2(|fft2(u1)|^2) = u1 ⊛ u1 = O⊛O + R*O(x+x0) + R*O*(x0-x) + R^2*delta。
% 思考题:撤掉回放透镜、眼睛放在合适位置只能看到一阶像的一部分(衍射阶
%     分离才行);整个过程就是"物 + 参考点"的自相关实现。
