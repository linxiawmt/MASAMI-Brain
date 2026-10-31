# ===== 书 2.5.2 · 离轴点源全息：记录（式 2.5-7）与再现（式 2.5-8） =====
# 改编自你已有的 holography.m（规范化重写 + 统一到书约定）：
# 点源位于 (0,0,z0)（胶片在 z=0，点源在 z=-z0），离轴平面参考波以记录角 theta 入射。
# 记录：t(x,y) = |psi_r + psi_o|^2 = A + B sin[ (k0/2z0)(x^2+y^2) + k0 x sin(theta) ]
# 再现：用同一参考波照射全息图，菲涅耳传播：
#     零级波沿 theta 方向继续出射                       (2.5-8a)
#     实像（孪生像）在胶片前方 z=+z0，横向 x=-2 z0 sin(theta)   (2.5-8b)
#     虚像（原始点源位置）在胶片后方 z=-z0，横向 (0,0)          (2.5-8c)
#
# MATLAB / Python 差异提示：
# 1) 本脚本全程严格书约定：前向变换 ifft2*N^2、逆变换 fft2/N^2，传播子用书式
#    (2.3-13) 的 exp(-j k0 z) exp(+j k^2 z / 2k0)。你旧版 holography.m 用 fft2 配
#    exp(+j k^2 z / 2k0)：fft2 与书约定互为镜像，两块强度图依然正确（相当于把
#    结果左右镜像），但复振幅是镜像的，且实/虚像的标注要与镜像约定配套才对——
#    本版用书协议后，实像严格落在 x = -2 z0 sin(theta)（式 2.5-8b），无歧义。
# 2) 物波按书 (2.5-7) 取 psi_o = (j k0 / 2 pi z0) exp[-j k0 r^2/(2 z0)]（省略全局
#    相位 exp(-j k0 z0)，书同样省略）；参考波振幅 a 取与 |psi_o| 最大值相同，
#    使条纹对比度 B/A 最大。
# 3) 采样检查：载频 fc = sin(theta)/lambda0，本例约 158 cycle/cm，dx≈9.8 um，
#    每条条纹约 6.5 个采样点，满足 Nyquist（见 2.5.3 数字全息的分辨率讨论）。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 物理参数 -----
L = 1.0                     # 胶片/探测器尺寸：1 cm x 1 cm
N = 1024                    # 采样点数（1024 便于解析载频条纹）
dx = L / (N - 1)
X = (np.arange(N)) * dx - L / 2
Y = (np.arange(N)) * dx - L / 2

lam = 0.6328e-4             # He-Ne 红光波长（cm）
k0 = 2 * np.pi / lam        # 波数
z0 = 10.0                   # 点源到胶片的距离（cm）
theta = 0.01                # 记录角（rad）；2*z0*sin(theta)=0.2 cm 在窗口内

# ----- 2. 频域网格（DFT 真实间距）与书式变换协议 -----
dk = 2 * np.pi / (N * dx)
kx = (np.arange(N) - N // 2) * dk
ky = (np.arange(N) - N // 2) * dk
KX, KY = np.meshgrid(kx, ky)

def ft2_book(U):
    """书式前向变换（+j 指数）= ifft2 * N^2"""
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    """书式逆变换（-j 指数）= fft2 / N^2"""
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

def prop_book(U0, z):
    """角谱传播（式 2.3-5 精确形式，自动滤消逝波）：z 可正可负"""
    root = 1 - (KX**2 + KY**2) / k0**2
    root[root < 0] = 0
    H = np.exp(-1j * k0 * z * np.sqrt(root))
    return ift2_book(ft2_book(U0) * H)

# ----- 3. 记录（Coding）：物波球面波 + 离轴参考平面波 -----
r2 = X[:, None]**2 + Y[None, :]**2
psi_o = (1j * k0 / (2 * np.pi * z0)) * np.exp(-1j * k0 * r2 / (2 * z0))   # 物波（发散球面波）
a = np.abs(psi_o).max()                                                  # 参考波振幅
psi_r = a * np.exp(1j * k0 * X * np.sin(theta))                          # 离轴参考波（书式 +j）

Hologram = np.abs(psi_r + psi_o)**2                                       # 全息图 t(x,y)
A0 = a**2 + (k0 / (2 * np.pi * z0))**2                                    # 式 (2.5-7) 的 A
B0 = 2 * a * k0 / (2 * np.pi * z0)                                        # 式 (2.5-7) 的 B = a k0/(pi z0)

# ----- 4. 再现（Decoding）：参考波照射 + 菲涅耳传播到 +z0 / -z0 -----
U_trans = Hologram * psi_r                       # psi_rc * t(x,y)
U_real = prop_book(U_trans, +z0)                 # 传播 +z0：实像（孪生像）平面
U_virtual = prop_book(U_trans, -z0)              # 反向传播 -z0：虚像平面
I_real = np.abs(U_real)**2
I_virtual = np.abs(U_virtual)**2

x_real_theory = -2 * z0 * np.sin(theta)          # 式 (2.5-8b)：delta(x + 2 z0 sinθ, y)
pk_r = np.unravel_index(np.argmax(I_real), I_real.shape)
pk_v = np.unravel_index(np.argmax(I_virtual), I_virtual.shape)
print(f'实像峰值位置  x = {X[pk_r[1]]:+.4f} cm （理论 {x_real_theory:+.4f} cm, y=0）')
print(f'虚像峰值位置  x = {X[pk_v[1]]:+.4f} cm （理论 {0:+.4f} cm, y=0）')

# ----- 5. 可视化 -----
fig, ax = plt.subplots(2, 2, figsize=(11.5, 9))
cx = np.abs(X) < 0.06                            # 全息图中心 1.2mm 放大看载频条纹

ax[0, 0].imshow(Hologram[np.ix_(cx, cx)], extent=[X[cx][0], X[cx][-1], Y[cx][0], Y[cx][-1]],
                cmap='gray', origin='lower')
ax[0, 0].set_title('Off-axis hologram t(x,y), zoomed center')
ax[0, 0].set_xlabel('x (cm)'); ax[0, 0].set_ylabel('y (cm)')

ax[0, 1].imshow(I_real**0.35, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[0, 1].set_title(f'Real image (twin) at z=+{z0:g} cm, x={X[pk_r[1]]:+.3f} cm')
ax[0, 1].set_xlabel('x (cm)'); ax[0, 1].set_ylabel('y (cm)')

ax[1, 0].imshow(I_virtual**0.35, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[1, 0].set_title(f'Virtual image at z=-{z0:g} cm, x={X[pk_v[1]]:+.3f} cm')
ax[1, 0].set_xlabel('x (cm)'); ax[1, 0].set_ylabel('y (cm)')

prof_r = I_real[N // 2, :] / I_real.max()
prof_v = I_virtual[N // 2, :] / I_virtual.max()
ax[1, 1].plot(X, prof_r, 'r-', label='real-image plane (z=+z0)')
ax[1, 1].plot(X, prof_v, 'b--', label='virtual-image plane (z=-z0)')
ax[1, 1].axvline(x_real_theory, color='k', ls=':', lw=1)
ax[1, 1].set_xlim(-0.4, 0.4); ax[1, 1].legend(); ax[1, 1].grid(True)
ax[1, 1].set_title('1-D intensity slices (y=0)')

plt.tight_layout()
plt.show()
