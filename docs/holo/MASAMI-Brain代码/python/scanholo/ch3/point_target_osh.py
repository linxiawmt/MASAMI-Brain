# ===== 书 3.5 节 · 点目标扫描全息：sin/cos FZP 全息图与复全息图重建（式 3.5-1 ~ 3.5-10） =====
# 物理设置：p1(x,y)=1（光瞳全开 -> L1 后焦面上点源，频率 w0）、p2(x,y)=delta（-> 平面波，频率 w0+Omega）。
# 点目标 |Gamma_0|^2 = delta(x,y)delta(z-z0) 位于距扫描镜 z0 处。
# 主结果：
#   H_sin(x,y) = (k0/2*pi*z0) sin[(k0/2*z0)(x^2+y^2)]           (3.5-4a)
#   H_cos(x,y) = (k0/2*pi*z0) cos[(k0/2*z0)(x^2+y^2)]           (3.5-4b)
#   H_c+ = H_cos + j*H_sin = j*i_Omegap  ->  与 h(x,y;z0) 卷积聚焦成点（conc = 1.0000）
#   H_c- = H_cos - j*H_sin = -j*conj(i_Omegap) -> 同样卷积得到 2*z0 处菲涅耳衍射图（corr = 1.0000）
#   sin 全息图单独重建：主像只占一半能量（另一半给孪生像/离焦斑，实测 0.501）——孪生像问题。
#
# 约定（全书）：前向变换 +j 指数 -> 用 ifft2*N^2；逆变换 -j 指数 -> 用 fft2/N^2。
# OTF_osh(kx,ky;z) = exp[-j*z/(2*k0)*(kx^2+ky^2)]  (3.5-1a)，|OTF|=1（纯相位）。
# 电子学符号：i_c = Re{i_Omegap}（cos 解调）、i_s = -Im{i_Omegap}（sin 解调，注意负号，
#             与式 (3.2-3b) 一致；原书 (3.4-14b) 直接写 Im，见笔记勘误）。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 1-based、Python 0-based：中心像素 MATLAB (257,257) <-> Python (256,256)。
# 2) imshow 需 origin='lower' 才与数学坐标一致；MATLAB imagesc 配 axis image。
# 3) e^{-j sigma k^2} 不是带限函数：离散记录的 sin/cos 全息图含 Poisson 旁瓣混叠
#    （书 fig 3.15b/c 颗粒状外观的数值根源），但闭环重建用 conj(OTFosh) 精确抵消，不受影响。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 物理参数（真实单位） -----
lam = 0.6328e-6            # He-Ne 红光波长（m）
k0 = 2*np.pi/lam           # 波数
z0 = 0.30                  # 点目标到扫描镜的距离（m）
N = 512                    # 网格
dx = 10e-6                 # 像素尺寸（m）
sigma = z0/(2*k0)          # = z0/(2*k0)，式 (3.5-1a) 的二次相位系数

x = (np.arange(N)-N//2)*dx
XX, YY = np.meshgrid(x, x)
df = 1/(N*dx)                                  # DFT 频率格距
kx = (np.arange(N)-N//2)*df
KY, KX = np.meshgrid(kx, kx)

# ----- 2. 书式傅里叶协议（前向 +j = ifft2*N^2，逆向 -j = fft2/N^2） -----
def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

# ----- 3. OSH 传递函数 (3.5-1a) 与点目标记录 -----
OTFosh = np.exp(-1j*sigma*(KX**2+KY**2))
I_pt = np.zeros((N, N)); I_pt[N//2, N//2] = 1.0
iOmegap = ift2_book(ft2_book(I_pt)*OTFosh)     # i_Omegap = F^{-1}{F{|Gamma0|^2} OTF_osh}  (3.4-15)

Hsin = np.real(iOmegap)                        # H_sin = Re{i_Omegap}   (3.5-3a)
Hcos = -np.imag(iOmegap)                       # H_cos = -Im{i_Omegap}  (3.5-3b，电子学负号)
Hc_plus  = Hcos + 1j*Hsin                      # H_c+ = H_cos + j*H_sin = j*i_Omegap    (3.5-5)
Hc_minus = Hcos - 1j*Hsin                      # H_c- = H_cos - j*H_sin = -j*conj(i_Omegap)

# ----- 4. 数字重建（式 3.5-10：与 h(x,y;z0) 卷积 = 频域乘 OTF_osh*） -----
rec_cp = ift2_book(ft2_book(Hc_plus)*np.conj(OTFosh))    # 聚焦 -> 点（实像）
rec_cm = ift2_book(ft2_book(Hc_minus)*np.conj(OTFosh))   # 离焦 -> 2*z0 菲涅耳衍射
rec_sp = ift2_book(ft2_book(Hsin)*np.conj(OTFosh))       # sin 全息图单独 -> 点 + 孪生

rcp = np.abs(rec_cp)**2
pk = np.unravel_index(np.argmax(rcp), rcp.shape)
conc = rcp[pk]/rcp.sum()
rsm = np.abs(rec_sp)**2
ratio_sin = rsm[pk]/rsm.sum()

# ----- 5. 对照：H_c- 重建 vs 原物在 2*z0 处的菲涅耳衍射（书 p.88 文字结论） -----
H2 = np.exp(-1j*k0*2*z0)*np.exp(1j*(KX**2+KY**2)*(2*z0)/(2*k0))
fres2 = ift2_book(ft2_book(I_pt)*H2)

def cc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return abs(np.vdot(a, b))/np.sqrt(np.vdot(a, a).real*np.vdot(b, b).real)

print(f'H_c+ 重建：峰值位置 {pk}（点目标在 {(N//2, N//2)}），能量集中度 = {conc:.4f}')
print(f'sin 全息图单独重建：主像能量占比 = {ratio_sin:.3f}（另一半进入孪生像）')
print(f'|H_c- 重建| vs |2*z0 菲涅耳衍射|：相关系数 = {cc(np.abs(rec_cm), np.abs(fres2)):.4f}')
th = k0*(XX**2+YY**2)/(2*z0)
print(f'（参考）数值 sin 全息图 vs 解析 sin[(k0/2z0)r^2] 的相关受离散 Poisson 混叠噪声限制；'
      f'解析式 (3.5-4a/b) 本身由 i_Omegap 的实/虚部直接给出')

# ----- 6. 可视化（对应笔记图 sim-point-target.png） -----
fig, ax = plt.subplots(2, 3, figsize=(13.5, 8))
ext = [x[0]*1e3, x[-1]*1e3, x[0]*1e3, x[-1]*1e3]
ax[0,0].imshow(np.sin(th), extent=ext, cmap='gray')
ax[0,0].set_title('sin-FZP of a point (Eq. 3.5-4a, analytic)')
ax[0,1].imshow(np.cos(th), extent=ext, cmap='gray')
ax[0,1].set_title('cos-FZP of a point (Eq. 3.5-4b, analytic)')
ax[0,2].imshow(rcp**0.3, extent=ext, cmap='gray')
ax[0,2].set_title(f'H_c+ * h(z0): focused point (conc = {conc:.2f})')
ax[1,0].imshow(rsm**0.3, extent=ext, cmap='gray')
ax[1,0].set_title('sin-FZP alone: point + twin (50/50)')
ax[1,1].imshow(np.abs(rec_cm)**0.5, extent=ext, cmap='gray')
ax[1,1].set_title('H_c- * h(z0): defocused')
ax[1,2].imshow(np.abs(fres2)**0.5, extent=ext, cmap='gray')
ax[1,2].set_title('Fresnel diffraction @2z0 (match)')
for a in ax.flat:
    a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
plt.tight_layout(); plt.show()
