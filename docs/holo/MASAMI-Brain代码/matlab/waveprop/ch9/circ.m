function out = circ(xa, ya, R)
% 书附录 B 的圆域函数（复现，ch9 依赖）
% circ(xa, ya) = 1，当 sqrt(xa^2 + ya^2) < 1/2（默认半径 1/2，直径 1）
% 三参数版 circ(xa, ya, R)：半径 R。Listing 9.8 里
% mask = circ(xn/D2, yn/D2, 1) 即"直径 D2 的望远镜口径"。
if nargin < 3
    R = 0.5;
end
out = double(sqrt(xa.^2 + ya.^2) <= R);
end
