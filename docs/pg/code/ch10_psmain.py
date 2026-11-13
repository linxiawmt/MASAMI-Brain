"""极客项目 ch10：粒子系统主程序（原书10章）"""
# 依赖：PyOpenGL + glfw + numpy + glutils.py（原书 common，见 ch09_simpleglfw.py 内附）；与 ch10_ps.py 同目录运行，需 star.png

import sys
import os
import math
from OpenGL.GL import *
import numpy
from ch10_ps import ParticleSystem, Camera  # 原书为 from ps import ...
import glutils
import glfw


class PSMaker:
    """GLFW Rendering window class for Particle System"""
    def __init__(self):
        self.camera = Camera([15.0, 0.0, 2.5],   # 眼睛：x=15、与喷泉同高
                             [0.0, 0.0, 2.5],    # 注视点：喷泉竖直轴上
                             [0.0, 0.0, 1.0])    # 向上向量：+z
        self.aspect = 1.0
        self.numP = 300
        self.t = 0  # 主程序时间（毫秒），只用于计时重启
        self.rotate = True

        cwd = os.getcwd()  # glfwInit 会改变当前工作目录，先保存
        glfw.glfwInit()
        os.chdir(cwd)

        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MAJOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MINOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_PROFILE, glfw.GLFW_OPENGL_CORE_PROFILE)

        self.width, self.height = 640, 480
        self.aspect = self.width/float(self.height)
        self.win = glfw.glfwCreateWindow(self.width, self.height, b"Particle System")
        glfw.glfwMakeContextCurrent(self.win)

        glViewport(0, 0, self.width, self.height)
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.2, 0.2, 0.2, 1.0)  # 深灰背景：加法混合的"底色"

        glfw.glfwSetMouseButtonCallback(self.win, self.onMouseButton)
        glfw.glfwSetKeyCallback(self.win, self.onKeyboard)
        glfw.glfwSetWindowSizeCallback(self.win, self.onSize)

        self.psys = ParticleSystem(self.numP)
        self.box = Box(1.0)
        self.exitNow = False

    def onMouseButton(self, win, button, action, mods):
        pass

    def onKeyboard(self, win, key, scancode, action, mods):
        if action == glfw.GLFW_PRESS:
            if key == glfw.GLFW_KEY_ESCAPE:  # ESC 退出
                self.exitNow = True
            elif key == glfw.GLFW_KEY_R:     # R：相机旋转开关
                self.rotate = not self.rotate
            elif key == glfw.GLFW_KEY_B:     # B：公告板开关
                self.psys.enableBillboard = not self.psys.enableBillboard
            elif key == glfw.GLFW_KEY_D:     # D：深度遮罩开关
                self.psys.disableDepthMask = not self.psys.disableDepthMask
            elif key == glfw.GLFW_KEY_T:     # T：透明混合开关
                self.psys.enableBlend = not self.psys.enableBlend

    def onSize(self, win, width, height):
        self.width, self.height = width, height
        self.aspect = width/float(height)
        glViewport(0, 0, self.width, self.height)

    def step(self):
        self.t += 10
        self.psys.step()
        if self.rotate:
            self.camera.rotate()
        if not int(self.t) % 5000:  # 每 5 秒重启（换随机颜色）
            self.psys.restart(self.numP)

    def run(self):
        glfw.glfwSetTime(0)
        t = 0.0
        while not glfw.glfwWindowShouldClose(self.win) and not self.exitNow:
            currT = glfw.glfwGetTime()
            if currT - t > 0.01:  # 每 0.01 秒渲染一帧
                t = currT
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                pMatrix = glutils.perspective(100.0, self.aspect, 0.1, 100.0)
                mvMatrix = glutils.lookAt(self.camera.eye, self.camera.center,
                                          self.camera.up)
                self.box.render(pMatrix, mvMatrix)  # 不透明物体先画
                self.psys.render(pMatrix, mvMatrix, self.camera)  # 半透明最后画
                self.step()
                glfw.glfwSwapBuffers(self.win)
            glfw.glfwPollEvents()
        glfw.glfwTerminate()


# ===== ch10 盒子 =====（原书 box.py：无光照暗红立方体，标出喷泉源头）
strVS = """
#version 330 core
in vec3 aVert;
uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;
out vec4 vCol;
void main() {
    gl_Position = uPMatrix * uMVMatrix * vec4(aVert, 1.0);
    vCol = vec4(0.8, 0.0, 0.0, 1.0);
}
"""

strFS = """
#version 330 core
in vec4 vCol;
out vec4 fragColor;
void main() {
    fragColor = vCol;
}
"""


class Box:
    def __init__(self, side):
        self.side = side

        self.program = glutils.loadShaders(strVS, strFS)
        glUseProgram(self.program)

        s = side/2.0
        vertices = [
            -s, s, -s,
            -s, -s, -s,
            s, s, -s,
            s, -s, -s,
            s, s, -s,
            -s, -s, -s,

            -s, s, s,
            -s, -s, s,
            s, s, s,
            s, -s, s,
            s, s, s,
            -s, -s, s,

            -s, -s, s,
            -s, -s, -s,
            s, -s, s,
            s, -s, -s,
            s, -s, s,
            -s, -s, -s,

            -s, s, s,
            -s, s, -s,
            s, s, s,
            s, s, -s,
            s, s, s,
            -s, s, -s,

            -s, -s, s,
            -s, -s, -s,
            -s, s, s,
            -s, s, -s,
            -s, s, s,
            -s, -s, -s,

            s, -s, s,
            s, -s, -s,
            s, s, s,
            s, s, -s,
            s, s, s,
            s, -s, -s
        ]

        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        vertexData = numpy.array(vertices, dtype=numpy.float32)
        self.vertexBuffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glBufferData(GL_ARRAY_BUFFER, 4*len(vertexData), vertexData, GL_STATIC_DRAW)
        self.vertIndex = glGetAttribLocation(self.program, b"aVert")
        glEnableVertexAttribArray(self.vertIndex)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glVertexAttribPointer(self.vertIndex, 3, GL_FLOAT, GL_FALSE, 0, None)
        glBindVertexArray(0)

    def render(self, pMatrix, mvMatrix):
        glUseProgram(self.program)
        glUniformMatrix4fv(glGetUniformLocation(self.program, b'uPMatrix'),
                           1, GL_FALSE, pMatrix)
        glUniformMatrix4fv(glGetUniformLocation(self.program, b'uMVMatrix'),
                           1, GL_FALSE, mvMatrix)
        glBindVertexArray(self.vao)
        glDrawArrays(GL_TRIANGLES, 0, 36)
        glBindVertexArray(0)


def main():
    print('starting particle system...')
    prog = PSMaker()
    prog.run()


if __name__ == '__main__':
    main()
