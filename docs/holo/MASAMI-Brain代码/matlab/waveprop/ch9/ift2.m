function g = ift2(G, delta_f)
% 书 Schmidt PM199 Listing 2.6：定标正确的二维逆 DFT（复现，ch9 依赖）
% g(x, y) = ∫∫ G(fx, fy) e^{+i2π(fx x + fy y)} dfx dfy 的离散近似：
% 乘 (N*delta_f)^2 对应 dfx dfy 的面积元。
g = fftshift(ifft2(fftshift(G))) * (size(G,1) * delta_f)^2;
end
