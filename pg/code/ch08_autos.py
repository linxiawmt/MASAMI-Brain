"""极客项目 ch08：三维立体画 Autostereogram（原书 8 章）"""

import sys, random, argparse
from PIL import Image, ImageDraw

def createRandomTile(dims):
    """create image filled with random dots"""
    img = Image.new('RGB', dims)
    draw = ImageDraw.Draw(img)
    r = int(min(*dims)/100)  # 圆点半径为短边的 1%
    n = 1000
    for i in range(n):
        # -r 让圆点保持在图内，平铺更干净
        x, y = random.randint(0, dims[0]-r), random.randint(0, dims[1]-r)
        fill = (random.randint(0, 255), random.randint(0, 255),
                random.randint(0, 255))
        draw.ellipse((x-r, y-r, x+r, y+r), fill)
    return img

def createTiledImage(tile, dims):
    """Create a larger image of size dims by tiling the given image"""
    img = Image.new('RGB', dims)
    W, H = dims
    w, h = tile.size
    cols = int(W/w) + 1
    rows = int(H/h) + 1
    for i in range(rows):
        for j in range(cols):
            img.paste(tile, (j*w, i*h))
    return img

def createAutostereogram(dmap, tile):
    """Given a depth map (image) and an input image, create a new image
    with pixels shifted according to depth"""
    if dmap.mode != 'L':  # Py3：is not 'L'→!=
        dmap = dmap.convert('L')
    # 未指定贴片则用随机圆点图
    if not tile:
        tile = createRandomTile((100, 100))
    # create an image by tiling
    img = createTiledImage(tile, dmap.size)
    sImg = img.copy()
    pixD = dmap.load()
    pixS = sImg.load()
    # shift pixels output based on depth map
    cols, rows = sImg.size
    for j in range(rows):
        for i in range(cols):
            xshift = pixD[i, j]//10  # Py3：/10 得浮点不能作像素索引
            xpos = i - tile.size[0] + xshift
            if xpos > 0 and xpos < cols:
                pixS[i, j] = pixS[xpos, j]
    return sImg

def main():
    print('creating autostereogram...')
    parser = argparse.ArgumentParser(description="Autosterograms...")
    parser.add_argument('--depth', dest='dmFile', required=True)
    parser.add_argument('--tile', dest='tileFile', required=False)
    parser.add_argument('--out', dest='outFile', required=False)
    args = parser.parse_args()

    outFile = 'as.png'
    if args.outFile:
        outFile = args.outFile

    tileFile = False
    if args.tileFile:
        tileFile = Image.open(args.tileFile)

    dmImg = Image.open(args.dmFile)
    asImg = createAutostereogram(dmImg, tileFile)
    asImg.save(outFile)
    print('done, wrote %s' % outFile)

if __name__ == '__main__':
    main()
