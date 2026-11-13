"""极客项目 ch09：理解 OpenGL——旋转贴图正方形、片段着色器剪圆（原书9章）"""
# 依赖：PyOpenGL + glfw + numpy + Pillow，OpenGL 3.3 core，无头环境跑不了；需同目录 star.png

import sys
import os
import math
from OpenGL.GL import *
from OpenGL.GL.shaders import compileShader
import numpy
from PIL import Image
import glfw


# ===== glutils =====（原书 common/glutils.py，本章用到四件，内联成单文件）
def loadShaders(strVS, strFS):
    """load vertex and fragment shaders from strings"""
    shaderV = compileShader([strVS], GL_VERTEX_SHADER)
    shaderF = compileShader([strFS], GL_FRAGMENT_SHADER)
    program = glCreateProgram()
    if not program:
        raise RuntimeError('glCreateProgram failed!')
    glAttachShader(program, shaderV)
    glAttachShader(program, shaderF)
    glLinkProgram(program)
    linked = glGetProgramiv(program, GL_LINK_STATUS)
    if not linked:
        glDeleteProgram(program)
        raise RuntimeError('Error linking program:\n%s\n' % glGetProgramInfoLog(program))
    return program


def perspective(fov, aspect, zNear, zFar):
    """returns matrix equivalent for gluPerspective"""
    fovR = math.radians(fov)
    f = 1.0/math.tan(fovR/2.0)
    return numpy.array([f/float(aspect), 0.0, 0.0, 0.0,
                        0.0, f, 0.0, 0.0,
                        0.0, 0.0, (zFar+zNear)/float(zNear-zFar), -1.0,
                        0.0, 0.0, 2.0*zFar*zNear/float(zNear-zFar), 0.0],
                       dtype=numpy.float32)  # 列优先铺平，与 OpenGL 一致


def lookAt(eye, center, up):
    """returns matrix equivalent of gluLookAt - based on MESA implementation"""
    m = numpy.identity(4, dtype=numpy.float32)
    forward = numpy.array(center) - numpy.array(eye)
    forward /= numpy.linalg.norm(forward)
    up /= numpy.linalg.norm(up)
    # side = forward x up；再重算 up = side x forward
    side = numpy.cross(forward, up)
    up = numpy.cross(side, forward)
    m[0][0] = side[0]
    m[1][0] = side[1]
    m[2][0] = side[2]
    m[0][1] = up[0]
    m[1][1] = up[1]
    m[2][1] = up[2]
    m[0][2] = -forward[0]
    m[1][2] = -forward[1]
    m[2][2] = -forward[2]
    # eye 平移
    t = numpy.identity(4, dtype=numpy.float32)
    t[3][0] += -eye[0]
    t[3][1] += -eye[1]
    t[3][2] += -eye[2]
    return t.dot(m)


def loadTexture(filename):
    """load OpenGL 2D texture from given image file"""
    img = Image.open(filename)
    imgData = numpy.array(list(img.getdata()), dtype=numpy.uint8)  # 原书 np.int8，NumPy 2.x 下溢出报 OverflowError，改 uint8
    texture = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, texture)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
    # 边缘截取 + 线性滤波（纹理坐标轴按惯例叫 S/T，不是 x/y）
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.size[0], img.size[1],
                 0, GL_RGBA, GL_UNSIGNED_BYTE, imgData)
    return texture


strVS = """
#version 330 core

layout(location = 0) in vec3 aVert;

uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;
uniform float uTheta;

out vec2 vTexCoord;

void main() {
    // rotational transform
    mat4 rot = mat4(
        vec4( cos(uTheta),  sin(uTheta), 0.0, 0.0),
        vec4(-sin(uTheta),  cos(uTheta), 0.0, 0.0),
        vec4( 0.0,          0.0,         1.0, 0.0),
        vec4( 0.0,          0.0,         0.0, 1.0)
    );
    // transform vertex
    gl_Position = uPMatrix * uMVMatrix * rot * vec4(aVert, 1.0);
    // set texture coord
    vTexCoord = aVert.xy + vec2(0.5, 0.5);
}
"""
strFS = """
#version 330 core

in vec2 vTexCoord;

uniform sampler2D tex2D;
uniform bool showCircle;

out vec4 fragColor;

void main() {
    if (showCircle) {
        // discard fragment outside circle
        if (distance(vTexCoord, vec2(0.5, 0.5)) > 0.5) {
            discard;
        }
        else {
            fragColor = texture(tex2D, vTexCoord);
        }
    }
    else {
        fragColor = texture(tex2D, vTexCoord);
    }
}
"""


class Scene:
    """OpenGL 3D scene class"""
    def __init__(self):
        # create shader
        self.program = loadShaders(strVS, strFS)
        glUseProgram(self.program)

        # 查着色器变量的位置
        self.pMatrixUniform = glGetUniformLocation(self.program, b'uPMatrix')
        self.mvMatrixUniform = glGetUniformLocation(self.program, b'uMVMatrix')
        self.tex2D = glGetUniformLocation(self.program, b'tex2D')

        # 三角形带顶点：以原点为中心、边长 1.0 的正方形
        vertexData = numpy.array(
            [-0.5, -0.5, 0.0,
             0.5, -0.5, 0.0,
             -0.5, 0.5, 0.0,
             0.5, 0.5, 0.0], dtype=numpy.float32)

        # set up vertex array object (VAO)
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vertexBuffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glBufferData(GL_ARRAY_BUFFER,
                     4*len(vertexData),  # 字节数：12 个 float × 4 字节
                     vertexData, GL_STATIC_DRAW)
        glEnableVertexAttribArray(0)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glBindVertexArray(0)

        self.t = 0  # 旋转角度（度）
        self.texId = loadTexture('star.png')
        self.showCircle = False

    def step(self):
        self.t = (self.t + 1) % 360
        glUniform1f(glGetUniformLocation(self.program, b'uTheta'),
                    math.radians(self.t))

    def render(self, pMatrix, mvMatrix):
        glUseProgram(self.program)
        glUniformMatrix4fv(self.pMatrixUniform, 1, GL_FALSE, pMatrix)
        glUniformMatrix4fv(self.mvMatrixUniform, 1, GL_FALSE, mvMatrix)
        glUniform1i(glGetUniformLocation(self.program, b'showCircle'),
                    self.showCircle)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self.texId)
        glUniform1i(self.tex2D, 0)
        glBindVertexArray(self.vao)
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4)
        glBindVertexArray(0)


class RenderWindow:
    """GLFW Rendering window class"""
    def __init__(self):
        cwd = os.getcwd()  # glfwInit 会改变当前工作目录，先保存
        glfw.glfwInit()
        os.chdir(cwd)

        # 版本 hints：OpenGL 3.3 core（必须在创建窗口之前设置）
        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MAJOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MINOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_PROFILE, glfw.GLFW_OPENGL_CORE_PROFILE)

        self.width, self.height = 640, 480
        self.aspect = self.width/float(self.height)
        self.win = glfw.glfwCreateWindow(self.width, self.height, b'simpleglfw')
        glfw.glfwMakeContextCurrent(self.win)

        glViewport(0, 0, self.width, self.height)
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.5, 0.5, 0.5, 1.0)

        glfw.glfwSetMouseButtonCallback(self.win, self.onMouseButton)
        glfw.glfwSetKeyCallback(self.win, self.onKeyboard)
        glfw.glfwSetWindowSizeCallback(self.win, self.onSize)

        self.scene = Scene()
        self.exitNow = False

    def onMouseButton(self, win, button, action, mods):
        pass

    def onKeyboard(self, win, key, scancode, action, mods):
        if action == glfw.GLFW_PRESS:
            if key == glfw.GLFW_KEY_ESCAPE:  # ESC 退出
                self.exitNow = True
            else:
                self.scene.showCircle = not self.scene.showCircle  # 其他键：翻转圆形剪裁

    def onSize(self, win, width, height):
        self.width, self.height = width, height
        self.aspect = width/float(height)
        glViewport(0, 0, self.width, self.height)

    def run(self):
        glfw.glfwSetTime(0)
        t = 0.0
        while not glfw.glfwWindowShouldClose(self.win) and not self.exitNow:
            currT = glfw.glfwGetTime()
            if currT - t > 0.1:  # 每 0.1 秒重绘一帧
                t = currT
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                pMatrix = perspective(45.0, self.aspect, 0.1, 100.0)
                mvMatrix = lookAt([0.0, 0.0, -2.0], [0.0, 0.0, 0.0],
                                  [0.0, 1.0, 0.0])
                self.scene.render(pMatrix, mvMatrix)
                self.scene.step()
                glfw.glfwSwapBuffers(self.win)
            glfw.glfwPollEvents()
        glfw.glfwTerminate()
        # （原书 9.5 清单在 run() 后另有一个无人调用的 step()：glfw.SwapBuffers 缺前缀，
        #   系排版残留的死代码，未收录）


def main():
    print("Starting simpleglfw. "
          "Press any key to toggle cut. Press ESC to quit.")
    rw = RenderWindow()
    rw.run()


if __name__ == '__main__':
    main()
