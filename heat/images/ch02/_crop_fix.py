# 修边 v2：命中选择 = y0 大于裁切框顶部的最近命中（即图注本身）
import fitz, os
p = r'E:\Desktop\Theodore L. Bergman, Adrienne S. Lavine, Frank P. Incropera, David P. DeWitt - Fundamentals of Heat and Mass Transfer, Seventh Edition  -John Wiley & Sons, Inc. (2011).pdf'
outdir = r'E:\Desktop\深度学习\mynotes\docs\heat\images\ch02'
doc = fitz.open(p)

def cap_below(pno, needle, ymin):
    hits = [h for h in doc[pno].search_for(needle) if h.y0 > ymin]
    assert hits, (pno, needle, ymin)
    return min(hits, key=lambda h: h.y0)

JOBS = [
 # name, pno, x0, y0, x1(None=keep), y1(None=keep), (needle for y1 cap), (needle,off for x1)
 ('fig-2-7',  99, 125, 495, None, 730, None, ('FIGURE 2.7', -4)),
 ('fig-2-10', 103, 190, 68, 470, None, ('FIGURE 2.10', -3), None),
 ('fig-2-11', 108, 168, 492, 500, None, ('FIGURE 2.11', -3), None),
 ('fig-2-12', 111, 170, 65, 495, None, ('FIGURE 2.12', -3), None),
 ('fig-2-13', 111, 168, 490, 495, None, ('FIGURE 2.13', -3), None),
]
for name, pno, x0, y0, x1, y1, yneedle, xadj in JOBS:
    if yneedle:
        y1 = cap_below(pno, yneedle, y0 + 30).y0 - 3
    if xadj:
        x1 = cap_below(pno, xadj[0], y0 + 30).x0 + xadj[1]
    page = doc[pno]
    r = fitz.Rect(x0, y0, x1, y1) & page.rect
    pix = page.get_pixmap(dpi=180, clip=r)
    pix.save(os.path.join(outdir, name + '.png'))
    print(name, [round(v,1) for v in r], pix.width, pix.height)

# e3: 顶部含 'A = 10 m²'，底部止于 Assumptions 前
b_top = doc[113].search_for('A = 10')[0]
b_bot = min(h for h in doc[113].search_for('One-dimensional conduction') if h.y0 > 300)
r = fitz.Rect(195, b_top.y0 - 6, 450, b_bot.y0 - 4)
doc[113].get_pixmap(dpi=180, clip=r).save(os.path.join(outdir, 'fig-2-e3.png'))
print('fig-2-e3', [round(v,1) for v in r])

# e4c / e4q / t22
b = min(h for h in doc[118].search_for('Note how the distributions') if h.y0 > 300)
r = fitz.Rect(298, 310, 545, b.y0 - 4)
doc[118].get_pixmap(dpi=180, clip=r).save(os.path.join(outdir, 'fig-2-e4c.png'))
print('fig-2-e4c', [round(v,1) for v in r])

b = min(h for h in doc[118].search_for('Convince yourself') if h.y0 > 500)
r = fitz.Rect(290, 512, 432, b.y0 - 4)
doc[118].get_pixmap(dpi=180, clip=r).save(os.path.join(outdir, 'fig-2-e4q.png'))
print('fig-2-e4q', [round(v,1) for v in r])

b = min(h for h in doc[116].search_for('diffusion equation at the surface') if h.y0 > 60)
r = fitz.Rect(335, b.y1 + 4, 430, 404)
doc[116].get_pixmap(dpi=180, clip=r).save(os.path.join(outdir, 'fig-2-t22.png'))
print('fig-2-t22', [round(v,1) for v in r])
