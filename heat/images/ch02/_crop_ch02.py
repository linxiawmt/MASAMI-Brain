# 临时副本：crop_book_figs.py 适配本书全大写图注 "FIGURE 2.x"
# 仅改 captions() 的编号正则为 ^FIGURE；用后即删
import re, sys, os
import fitz

BOOK = r'E:\Desktop\Theodore L. Bergman, Adrienne S. Lavine, Frank P. Incropera, David P. DeWitt - Fundamentals of Heat and Mass Transfer, Seventh Edition  -John Wiley & Sons, Inc. (2011).pdf'
p0, p1, outdir, prefix = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
os.makedirs(outdir, exist_ok=True)
doc = fitz.open(BOOK)
DPI = 180

def captions(page):
    """返回 [(编号, caption行bbox)]"""
    out = []
    d = page.get_text('dict')
    for b in d['blocks']:
        for l in b.get('lines', []):
            txt = ''.join(s['text'] for s in l['spans']).strip()
            m = re.match(r'^FIGURE\s+(\d+\.\d+)\b', txt)
            if m:
                out.append((m.group(1), fitz.Rect(l['bbox'])))
    return out

def crop_for(page, cap):
    """给定图注 bbox，找其上方的图形区域并合并；无则 None"""
    cands = [fitz.Rect(c['rect']) if 'rect' in c else fitz.Rect(c['bbox'])
             for c in [dict(type='img', bbox=i['bbox']) for i in page.get_image_info()]]
    cands += [fitz.Rect(r) for r in page.cluster_drawings()]
    cands = [r for r in cands if r.width > 50 and r.height > 25]
    pool = [r for r in cands
            if r.y1 <= cap.y0 + 10 and (cap.y0 - r.y1) < 300
            and min(r.x1, cap.x1) - max(r.x0, cap.x0) > -20]
    if not pool:
        return None
    seed = max(pool, key=lambda r: r.y1)
    merged = fitz.Rect(seed)
    for _ in range(4):
        grow = fitz.Rect(merged.x0-25, merged.y0-25, merged.x1+25, merged.y1+25)
        for r in pool:
            if r.intersects(grow):
                merged |= r
    return merged

saved, missed = [], []
for pno in range(p0, p1+1):
    page = doc[pno]
    for num, cap in captions(page):
        rect = crop_for(page, cap)
        if rect is None:
            missed.append((num, pno+1)); continue
        rect = fitz.Rect(max(rect.x0-4,0), max(rect.y0-4,0),
                         min(rect.x1+4, page.rect.x1), min(rect.y1+3, page.rect.y1))
        pix = page.get_pixmap(dpi=DPI, clip=rect)
        fp = os.path.join(outdir, f'fig-{num.replace(".", "-")}.png')
        pix.save(fp); saved.append((num, pno+1))

print(f'saved {len(saved)}:', ' '.join(n for n,_ in saved))
if missed:
    print('MISSED:', missed)

from PIL import Image, ImageDraw
files = sorted(f for f in os.listdir(outdir) if f.startswith('fig-') and f.endswith('.png'))
if files:
    TH = 240
    cols = 4
    rows = (len(files)+cols-1)//cols
    thumbs = []
    for f in files:
        im = Image.open(os.path.join(outdir, f)); w,h = im.size
        im = im.resize((int(w*TH/h), TH))
        cv = Image.new('RGB', (im.width, TH+18), 'white')
        cv.paste(im, (0,0))
        d = ImageDraw.Draw(cv); d.text((4, TH+2), f, fill='black')
        thumbs.append(cv)
    W = max(t.width for t in thumbs)
    sheet = Image.new('RGB', (cols*(W+8), rows*(TH+22)), 'white')
    for i, t in enumerate(thumbs):
        sheet.paste(t, ((i%cols)*(W+8), (i//cols)*(TH+22)))
    sheet.save(os.path.join(outdir, f'_{prefix}_sheet.png'))
    print('sheet:', os.path.join(outdir, f'_{prefix}_sheet.png'))
