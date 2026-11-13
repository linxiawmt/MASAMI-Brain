"""极客项目 ch01：iTunes 播放列表分析（原书 1 章）"""

import os, argparse
import tempfile
from matplotlib import pyplot
import plistlib
import numpy as np

# 演示数据（曲名, 时长ms, 专辑评分）：含同名同时长重复、同名不同时长变体
DEMO_TRACKS = [
    ('Neon Sky', 245000, 80), ('Paper Moon', 199000, 60), ('Starlight', 231000, 100),
    ('Iron Rain', 254000, 40), ('Silver Birch', 187000, 100), ('Low Tide', 262000, 80),
    ('Iron Rain', 311000, 60), ('Neon Sky', 245000, 80), ('Glass River', 208000, 40),
    ('Cold Coffee', 293000, 20),
]

_demos = None

def makeDemoPlists():
    """生成的演示 plist 写入系统临时目录：demo1 全 10 音轨 2 播放列表，demo2 取 4-10 号音轨供 --common 交集。"""
    global _demos
    if _demos is None:
        d = tempfile.mkdtemp(prefix='pg_demo_')
        tracks = {str(i): {'Track ID': i, 'Name': n, 'Total Time': t, 'Album Rating': r}
                  for i, (n, t, r) in enumerate(DEMO_TRACKS, 1)}
        plists = [{'Name': 'Mix A', 'Playlist Items': [{'Track ID': i} for i in range(1, 7)]},
                  {'Name': 'Mix B', 'Playlist Items': [{'Track ID': i} for i in range(5, 11)]}]
        top = {'Major Version': 1, 'Application Version': 'demo',
               'Tracks': tracks, 'Playlists': plists}
        f1 = os.path.join(d, 'demo1.xml')
        plistlib.dump(top, open(f1, 'wb'))
        top2 = dict(top, Tracks={k: v for k, v in tracks.items() if int(k) >= 4}, Playlists=plists[1:])
        f2 = os.path.join(d, 'demo2.xml')
        plistlib.dump(top2, open(f2, 'wb'))
        print('-- 指定的 plist 不存在，已生成演示数据：%s' % d)
        _demos = [f1, f2]
    return _demos

def findCommonTracks(fileNames):
    """Find common tracks in given playlist files, and save them to common.txt."""
    trackNameSets = []
    for fileName in fileNames:
        trackNames = set()
        plist = plistlib.load(open(fileName, 'rb'))  # Py3：readPlist→load
        tracks = plist['Tracks']
        for trackId, track in tracks.items():
            try:
                trackNames.add(track['Name'])
            except KeyError:
                pass
        trackNameSets.append(trackNames)
    commonTracks = set.intersection(*trackNameSets)
    if len(commonTracks) > 0:
        with open('common.txt', 'w', encoding='UTF-8') as f:  # Py3：文本模式直接写 str
            for val in commonTracks:
                f.write('%s\n' % val)
        print("%d common tracks found. Track names written to common.txt." % len(commonTracks))
    else:
        print("No common tracks!")

def plotStats(fileName):
    """Plot some statistics by reading track information from playlist."""
    plist = plistlib.load(open(fileName, 'rb'))  # Py3：readPlist→load
    tracks = plist['Tracks']
    ratings = []
    durations = []
    for trackId, track in tracks.items():
        try:
            ratings.append(track['Album Rating'])
            durations.append(track['Total Time'])
        except KeyError:
            pass
    if ratings == [] or durations == []:
        print("No valid Album Rating/Total Time data in %s." % fileName)
        return
    x = np.array(durations, np.int32)
    x = x/60000.0
    y = np.array(ratings, np.int32)
    pyplot.subplot(2, 1, 1)
    pyplot.plot(x, y, 'o')
    pyplot.axis([0, 1.05*np.max(x), -1, 110])
    pyplot.xlabel('Track duration')
    pyplot.ylabel('Track rating')
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x, bins=20)
    pyplot.xlabel('Track duration')
    pyplot.ylabel('Count')
    pyplot.show()

def findDuplicates(fileName):
    """Find duplicate tracks in given playlist."""
    print('Finding duplicate tracks in %s...' % fileName)
    plist = plistlib.load(open(fileName, 'rb'))  # Py3：readPlist→load
    tracks = plist['Tracks']
    trackNames = {}
    for trackId, track in tracks.items():
        try:
            name = track['Name']
            duration = track['Total Time']
            if name in trackNames:
                # 同名且时长（取整到秒）相同才计为重复
                if duration//1000 == trackNames[name][0]//1000:
                    count = trackNames[name][1]
                    trackNames[name] = (duration, count+1)
            else:
                trackNames[name] = (duration, 1)
        except KeyError:
            pass
    dups = []
    for k, v in trackNames.items():
        if v[1] > 1:
            dups.append((v[1], k))
    if len(dups) > 0:
        print("Found %d duplicates. Track names saved to dup.txt" % len(dups))
    else:
        print("No duplicate tracks found!")
    with open('dups.txt', 'w', encoding='UTF-8') as f:  # Py3：文本模式直接写 str
        for val in dups:
            f.write("[%d] %s\n" % (val[0], val[1]))

def main():
    descStr = "This program analyzes playlist files (.xml) exported from iTunes."
    parser = argparse.ArgumentParser(description=descStr)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--file', dest='plFile', required=False,
                       help='单个 plist：查找重复曲目并绘制统计图')
    group.add_argument('--common', nargs='*', dest='plFiles', required=False)
    group.add_argument('--stats', dest='plStats', required=False)
    group.add_argument('--dup', dest='plDup', required=False)
    args = parser.parse_args()

    # 文件不存在时改用自动生成的演示 plist，保证零依赖可跑
    if args.plFile and not os.path.exists(args.plFile):
        args.plFile = makeDemoPlists()[0]
    if args.plFiles:
        demos = makeDemoPlists()
        args.plFiles = [f if os.path.exists(f) else demos[i % 2]
                        for i, f in enumerate(args.plFiles)]
    if args.plStats and not os.path.exists(args.plStats):
        args.plStats = makeDemoPlists()[0]
    if args.plDup and not os.path.exists(args.plDup):
        args.plDup = makeDemoPlists()[0]

    if args.plFile:
        findDuplicates(args.plFile)
        plotStats(args.plFile)
    elif args.plFiles:
        findCommonTracks(args.plFiles)
    elif args.plStats:
        plotStats(args.plStats)
    elif args.plDup:
        findDuplicates(args.plDup)
    else:
        print("These are not the tracks you are looking for.")

if __name__ == '__main__':
    main()
