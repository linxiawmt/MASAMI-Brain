# 第2章最终裁图：手工定位框（pt），180dpi；另渲染公式核对区
import fitz, os
p = r'E:\Desktop\Theodore L. Bergman, Adrienne S. Lavine, Frank P. Incropera, David P. DeWitt - Fundamentals of Heat and Mass Transfer, Seventh Edition  -John Wiley & Sons, Inc. (2011).pdf'
outdir = r'E:\Desktop\深度学习\mynotes\docs\heat\images\ch02'
doc = fitz.open(p)

CROPS = {
 'fig-2-2':  (94, (138, 540, 252, 658)),
 'fig-2-3':  (95, (112, 74, 232, 158)),
 'fig-2-4':  (96, (125, 465, 380, 622)),
 'fig-2-5':  (97, (125, 358, 352, 648)),
 'fig-2-6':  (98, (138, 450, 357, 612)),
 'fig-2-7':  (99, (125, 495, 358, 732)),
 'fig-2-8':  (100, (140, 398, 372, 652)),
 'fig-2-9':  (101, (128, 396, 350, 645)),
 'fig-2-10': (103, (190, 68, 470, 250)),
 'fig-2-11': (108, (168, 492, 500, 702)),
 'fig-2-12': (111, (170, 65, 495, 268)),
 'fig-2-13': (111, (168, 490, 495, 665)),
 'fig-2-e2': (106, (255, 95, 460, 200)),
 'fig-2-e3': (113, (200, 295, 445, 485)),
 'fig-2-t22': (116, (335, 90, 430, 404)),
 'fig-2-e4': (117, (205, 110, 460, 280)),
 'fig-2-e4c': (118, (298, 310, 545, 440)),
 'fig-2-e4q': (118, (290, 512, 432, 648)),
}
for name, (pno, (x0, y0, x1, y1)) in CROPS.items():
    page = doc[pno]
    r = fitz.Rect(x0, y0, x1, y1) & page.rect
    pix = page.get_pixmap(dpi=180, clip=r)
    pix.save(os.path.join(outdir, name + '.png'))
    print(name, pix.width, pix.height)

# 公式核对区：搜索定位
for label, pno, needle in [('eq-2-12', 101, '(2.12)'), ('eq-2-30', 115, '(2.30)'),
                            ('eq-2-26', 112, '(2.26)'), ('eq-2-29', 112, '(2.29)'),
                            ('eq-2-7', 96, '(2.7)')]:
    page = doc[pno]
    hits = page.search_for(needle)
    if not hits:
        print('NOT FOUND', label); continue
    h = hits[-1]
    r = fitz.Rect(100, max(h.y0 - 42, 0), 512, min(h.y1 + 16, page.rect.y1))
    pix = page.get_pixmap(dpi=170, clip=r)
    fp = os.path.join(outdir, f'_eq_{label}.png')
    pix.save(fp)
    print(fp)
