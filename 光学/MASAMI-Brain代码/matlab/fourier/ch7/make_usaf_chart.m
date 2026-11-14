% ===== 书 7.2.3 节:make_usaf_chart 生成 USAF 1951 分辨率板 PNG =====
% 与 python/fourier/ch7/ch7_usaf_chart.py 同规格,替代原书网站已失效的
% USAF1951B250.png(书说可从 http://www.ece.nmsu.edu/~davvoelz/cfo/ 下载,
% 该地址现已不可用;如果你手里有原书 PNG,把它改名为 USAF1951B250.png
% 放进当前目录,coh_image.m 第 3 行的 imread 可原样使用)。
% 规格:250x250、黑底白条、映射 L=0.3 mm(du=1.2 um);
% USAF 频率标尺 f(g,n) = f00*2^(g+n/6),锚定 f00(组0/单元1)= 2e5 cyc/m,
% 使"组-2/单元6"= 1e5 cyc/m(周期 10 um),与书 7.2.3 节判读一致。
% 每单元 3 条竖条:条宽 1/(2f),条长 5 倍条宽;4 倍超采样抗混叠后二值化,
% 超出奈奎斯特频率的细条纹自然消失(不给理想图像混入假频率)。
% MATLAB/Python 差异:MATLAB 阵列 (row,col)、imwrite 直接写 uint8;
% Python 用 PIL 的 BOX 降采样,这里用块平均手写。

function make_usaf_chart()
Npx = 250;  Limg = 0.3e-3;  du = Limg/Npx;
F00 = 2.0e5;  SS = 4;
Ssz = Npx*SS;                 %注意:别用 size 当变量名(遮蔽内置函数)
img = zeros(Ssz, Ssz);
groups = [-2 -1 0 1 2];  elems = 1:6;
x_left = 30*SS;  y_top = 18*SS;
col_w = (Npx-30)/5*SS;

for i = 1:5
    g = groups(i);  x0 = x_left + (i-1)*col_w;  y = y_top;
    for n = elems
        f = F00 * 2^(g + n/6);       % 该单元空间频率 [cyc/m]
        w = 1/(2*f)/du*SS;           % 条宽(超采样像素)
        ln = 5*w;                    % 条长 = 5 倍条宽
        if y + ln > Ssz, break; end
        for b = 0:2                  % 3 条,间隔 = 条宽(周期 2w)
            bx = x0 + b*2*w;
            img(round(y):round(y+ln), round(bx):round(bx+w)) = 1;
        end
        y = y + ln + 4*SS;
    end
end

% 块平均降采样(等价抗混叠)+ 二值化
small = zeros(Npx, Npx);
for r = 1:Npx
    for c = 1:Npx
        blk = img((r-1)*SS+1 : r*SS, (c-1)*SS+1 : c*SS);
        small(r,c) = mean(blk(:));
    end
end
small = double(small > 0.5);
imwrite(uint8(255*small), 'USAF1951_py250.png');
fprintf('已生成 USAF1951_py250.png(250x250 二值)\n');
end
