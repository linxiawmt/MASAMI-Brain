"""极客项目 ch12：Arduino 光敏电阻（LDR）实时绘图（原书12章）"""
# 依赖：pyserial + matplotlib；python ch12_ldr.py --port COM3 [--N 100]（选项须在脚本名之后）

import serial
import argparse
from collections import deque

import matplotlib.pyplot as plt
import matplotlib.animation as animation


# plot class
class AnalogPlot:
    def __init__(self, strPort, maxLen):
        self.ser = serial.Serial(strPort, 9600)  # 波特率与 Arduino 端 Serial.begin 一致
        self.a0Vals = deque([0.0]*maxLen)
        self.a1Vals = deque([0.0]*maxLen)
        self.maxLen = maxLen

    # add data
    def add(self, data):
        assert len(data) == 2
        self.addToDeq(self.a0Vals, data[0])
        self.addToDeq(self.a1Vals, data[1])

    # add to deque, pop oldest value
    def addToDeq(self, buf, val):
        buf.pop()            # 从尾部删最老的
        buf.appendleft(val)  # 新值进头部 → 图形向右滚动

    # update plot
    def update(self, frameNum, a0, a1):
        try:
            line = self.ser.readline()
            data = [float(v) for v in line.decode().split()]  # Python 3：readline() 返回 bytes
            if len(data) == 2:
                self.add(data)
                a0.set_data(range(self.maxLen), self.a0Vals)
                a1.set_data(range(self.maxLen), self.a1Vals)
        except:
            pass  # 坏数据直接忽略
        return a0, a1

    # clean up
    def close(self):
        self.ser.flush()
        self.ser.close()


# main() function
def main():
    parser = argparse.ArgumentParser(description="LDR serial")
    parser.add_argument('--port', dest='port', required=True)
    parser.add_argument('--N', dest='maxLen', required=False)
    args = parser.parse_args()

    strPort = args.port
    print('reading from serial port %s...' % strPort)

    maxLen = 100
    if args.maxLen:
        maxLen = int(args.maxLen)

    analogPlot = AnalogPlot(strPort, maxLen)
    print('plotting data...')

    fig = plt.figure()
    ax = plt.axes(xlim=(0, maxLen), ylim=(0, 1023))
    a0, = ax.plot([], [])
    a1, = ax.plot([], [])
    anim = animation.FuncAnimation(fig, analogPlot.update,  # 存入 anim 变量：防回调被垃圾回收
                                   fargs=(a0, a1),
                                   interval=20)

    plt.show()
    analogPlot.close()
    print('exiting.')


if __name__ == '__main__':
    main()
