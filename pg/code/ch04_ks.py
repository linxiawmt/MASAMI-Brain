"""极客项目 ch04：Karplus-Strong 算法生成音乐（原书 4 章）"""

import sys, os
import time, random
import wave, argparse
import numpy as np
from collections import deque
from matplotlib import pyplot as plt

try:
    import pygame
except ImportError:
    pygame = None  # pip install pygame 后可用 --display/--play/--piano

# show plot of algorithm in action?
gShowPlot = False

# 五声音阶 C4-Eb-F-G-Bb
pmNotes = {'C4': 262, 'Eb': 311, 'F': 349, 'G': 391, 'Bb': 466}

def writeWAVE(fname, data):
    file = wave.open(fname, 'wb')
    nChannels = 1
    sampleWidth = 2
    frameRate = 44100
    nFrames = 44100
    file.setparams((nChannels, sampleWidth, frameRate, nFrames,
                    'NONE', 'noncompressed'))
    file.writeframes(data)
    file.close()

def generateNote(freq):
    nSamples = 44100
    sampleRate = 44100
    N = int(sampleRate/freq)  # 环形缓冲长度 = 一个周期
    buf = deque([random.random() - 0.5 for i in range(N)])
    if gShowPlot:
        axline, = plt.plot(buf)
    samples = np.array([0]*nSamples, 'float32')
    for i in range(nSamples):
        samples[i] = buf[0]
        avg = 0.995*0.5*(buf[0] + buf[1])  # 衰减平均
        buf.append(avg)
        buf.popleft()
        if gShowPlot:
            if i % 1000 == 0:
                axline.set_ydata(buf)
                plt.draw()
    # 16-bit PCM
    samples = np.array(samples * 32767, 'int16')
    return samples.tobytes()  # Py3：tostring()→tobytes()

# play a wav file
class NotePlayer:
    def __init__(self):
        pygame.mixer.pre_init(44100, -16, 1, 2048)
        pygame.init()
        self.notes = {}

    def add(self, fileName):
        self.notes[fileName] = pygame.mixer.Sound(fileName)

    def play(self, fileName):
        try:
            self.notes[fileName].play()
        except KeyError:
            print(fileName + ' not found!')

    def playRandom(self):
        """play a random note"""
        index = random.randint(0, len(self.notes)-1)
        note = list(self.notes.values())[index]
        note.play()

def main():
    global gShowPlot

    parser = argparse.ArgumentParser(description="Generating sounds with Karplus String Algorithm.")
    parser.add_argument('--display', action='store_true', required=False)
    parser.add_argument('--play', action='store_true', required=False)
    parser.add_argument('--piano', action='store_true', required=False)
    args = parser.parse_args()

    if (args.display or args.play or args.piano) and pygame is None:
        sys.exit('需要 pygame：pip install pygame 后可用 --display/--play/--piano')

    # show plot if flag set
    if args.display:
        gShowPlot = True
        plt.ion()

    nplayer = NotePlayer() if pygame is not None else None

    print('creating notes...')
    for name, freq in list(pmNotes.items()):
        fileName = name + '.wav'
        if not os.path.exists(fileName) or args.display:
            data = generateNote(freq)
            print('creating ' + fileName + '...')
            writeWAVE(fileName, data)
        else:
            print('fileName already created. skipping...')

        if nplayer is not None:
            nplayer.add(name + '.wav')

        if args.display:
            nplayer.play(name + '.wav')
            time.sleep(0.5)

    # play a random tune
    if args.play:
        while True:
            try:
                nplayer.playRandom()
                # rest - 1 to 8 beats
                rest = np.random.choice([1, 2, 4, 8], 1,
                                        p=[0.15, 0.7, 0.1, 0.05])
                time.sleep(0.25*rest[0])
            except KeyboardInterrupt:
                exit()

    # random piano mode
    if args.piano:
        while True:
            for event in pygame.event.get():
                if (event.type == pygame.KEYUP):
                    print("key pressed")
                    nplayer.playRandom()
                    time.sleep(0.5)

if __name__ == '__main__':
    main()
