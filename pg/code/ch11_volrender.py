"""极客项目 ch11：体渲染——体光线投射 + 三维切片（原书11章）"""
# 依赖：PyOpenGL + glfw + numpy + Pillow，OpenGL 3.3 core，无头环境跑不了；运行：python ch11_volrender.py --dir 切片图片目录

import sys
import os
import math
import argparse
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
                       dtype=numpy.float32)


def ortho(l, r, b, t, n, f):
    """returns matrix equivalent of glOrtho"""
    return numpy.array([2.0/float(r-l), 0.0, 0.0, 0.0,
                        0.0, 2.0/float(t-b), 0.0, 0.0,
                        0.0, 0.0, -2.0/float(f-n), 0.0,
                        -(r+l)/float(r-l), -(t+b)/float(t-b),
                        -(f+n)/float(f-n), 1.0],
                       dtype=numpy.float32)


def lookAt(eye, center, up):
    """returns matrix equivalent of gluLookAt - based on MESA implementation"""
    m = numpy.identity(4, dtype=numpy.float32)
    forward = numpy.array(center) - numpy.array(eye)
    forward /= numpy.linalg.norm(forward)
    up /= numpy.linalg.norm(up)
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
    t = numpy.identity(4, dtype=numpy.float32)
    t[3][0] += -eye[0]
    t[3][1] += -eye[1]
    t[3][2] += -eye[2]
    return t.dot(m)


# ===== volreader =====
def loadVolume(dirName):
    """read volume from directory as a 3D texture"""
    files = sorted(os.listdir(dirName))  # 必须排序：切片顺序即深度方向
    print('loading images from: %s' % dirName)
    imgDataList = []
    count = 0
    width, height = 0, 0
    for file in files:
        file_path = os.path.abspath(os.path.join(dirName, file))
        try:
            img = Image.open(file_path)
            imgData = numpy.array(img.getdata(), dtype=numpy.uint8)
            # 尺寸校验：全部切片必须一致
            if count == 0:
                width, height = img.size[0], img.size[1]
                imgDataList.append(imgData)
            elif (width, height) == (img.size[0], img.size[1]):
                imgDataList.append(imgData)
            else:
                print('mismatch')
                raise RuntimeError("image size mismatch")  # 原书拼作 RunTimeError，且被裸 except 吞掉
            count += 1
        except OSError:
            print('Invalid image: %s' % file_path)  # 非图像文件：跳过

    depth = count
    data = numpy.concatenate(imgDataList)
    print('volume data dims: %d %d %d' % (width, height, depth))

    # 上传为三维纹理（单通道 8 位，GL_RED）
    texture = glGenTextures(1)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)  # 灰度行宽不必 4 字节对齐
    glBindTexture(GL_TEXTURE_3D, texture)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RED,
                 width, height, depth, 0,
                 GL_RED, GL_UNSIGNED_BYTE, data)
    return (texture, width, height, depth)


# ===== raycube =====（原书 raycube.py；其 strVS/strFS 与另两模块重名，合并后按节改名）
cubeVS = """
#version 330 core

layout(location = 1) in vec3 cubePos;
layout(location = 2) in vec3 cubeCol;

uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;
out vec4 vColor;

void main()
{
    // set back face color
    vColor = vec4(cubeCol.rgb, 1.0);
    // transformed position
    vec4 newPos = vec4(cubePos.xyz, 1.0);
    // set position
    gl_Position = uPMatrix * uMVMatrix * newPos;
}
"""
cubeFS = """
#version 330 core

in vec4 vColor;
out vec4 fragColor;

void main()
{
    fragColor = vColor;
}
"""


class RayCube:
    """class used to generate rays used in ray casting"""

    def __init__(self, width, height):
        self.width, self.height = width, height
        self.program = loadShaders(cubeVS, cubeFS)

        # 立方体顶点：颜色 = 坐标（彩色立方体技巧）
        vertices = numpy.array([
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0
        ], dtype=numpy.float32)
        colors = numpy.array([
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0
        ], dtype=numpy.float32)

        # 12 个三角形（6 个面）的索引
        indices = numpy.array([
            4, 5, 7,
            7, 5, 6,
            5, 1, 6,
            6, 1, 2,
            1, 0, 2,
            2, 0, 3,
            0, 4, 3,
            3, 4, 7,
            6, 2, 7,
            7, 2, 3,
            4, 0, 5,
            5, 0, 1
        ], dtype=numpy.int16)

        self.nIndices = indices.size

        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)

        self.vertexBuffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glBufferData(GL_ARRAY_BUFFER, 4*len(vertices), vertices, GL_STATIC_DRAW)

        self.colorBuffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.colorBuffer)
        glBufferData(GL_ARRAY_BUFFER, 4*len(colors), colors, GL_STATIC_DRAW)

        # 索引缓冲（int16 → GL_UNSIGNED_SHORT）
        self.indexBuffer = glGenBuffers(1)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.indexBuffer)
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2*len(indices), indices, GL_STATIC_DRAW)

        # 属性位置与顶点着色器 layout(location) 对应
        aPosLoc = 1
        aColorLoc = 2
        glEnableVertexAttribArray(1)
        glEnableVertexAttribArray(2)

        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glVertexAttribPointer(aPosLoc, 3, GL_FLOAT, GL_FALSE, 0, None)

        glBindBuffer(GL_ARRAY_BUFFER, self.colorBuffer)
        glVertexAttribPointer(aColorLoc, 3, GL_FLOAT, GL_FALSE, 0, None)

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.indexBuffer)
        glBindVertexArray(0)

        self.initFBO()

    def renderBackFace(self, pMatrix, mvMatrix):
        """renders back-face of ray-cube to a texture and returns it"""
        glBindFramebuffer(GL_FRAMEBUFFER, self.fboHandle)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self.texHandle)
        # 剔除正面 → 只画背面
        self.renderCube(pMatrix, mvMatrix, self.program, True)
        glBindTexture(GL_TEXTURE_2D, 0)
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        glBindRenderbuffer(GL_RENDERBUFFER, 0)
        return self.texHandle

    def renderFrontFace(self, pMatrix, mvMatrix, program):
        """render front face of ray-cube"""
        self.renderCube(pMatrix, mvMatrix, program, False)

    def renderCube(self, pMatrix, mvMatrix, program, cullFace):
        """render cube use face culling if flag set"""
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glUseProgram(program)
        glUniformMatrix4fv(glGetUniformLocation(program, b'uPMatrix'),
                           1, GL_FALSE, pMatrix)
        glUniformMatrix4fv(glGetUniformLocation(program, b'uMVMatrix'),
                           1, GL_FALSE, mvMatrix)

        # cullFace=True：剔除正面、只渲染背面到 FBO
        glDisable(GL_CULL_FACE)
        if cullFace:
            glFrontFace(GL_CCW)
            glCullFace(GL_FRONT)
            glEnable(GL_CULL_FACE)

        glBindVertexArray(self.vao)
        glDrawElements(GL_TRIANGLES, self.nIndices, GL_UNSIGNED_SHORT, None)
        glBindVertexArray(0)

        if cullFace:
            glDisable(GL_CULL_FACE)

    def reshape(self, width, height):
        self.width = width
        self.height = height
        self.aspect = width/float(height)
        self.clearFBO()
        self.initFBO()

    def initFBO(self):
        # 帧缓冲对象：背面渲染结果存进纹理
        self.fboHandle = glGenFramebuffers(1)
        self.texHandle = glGenTextures(1)
        self.depthHandle = glGenRenderbuffers(1)

        glBindFramebuffer(GL_FRAMEBUFFER, self.fboHandle)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self.texHandle)

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, self.width, self.height,
                     0, GL_RGBA, GL_UNSIGNED_BYTE, None)

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D, self.texHandle, 0)

        glBindRenderbuffer(GL_RENDERBUFFER, self.depthHandle)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24,
                              self.width, self.height)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                  GL_RENDERBUFFER, self.depthHandle)

        status = glCheckFramebufferStatus(GL_FRAMEBUFFER)
        if status == GL_FRAMEBUFFER_COMPLETE:
            pass
        elif status == GL_FRAMEBUFFER_UNSUPPORTED:
            print("fbo %d unsupported" % self.fboHandle)
        else:
            print("fbo %d Error" % self.fboHandle)

        glBindTexture(GL_TEXTURE_2D, 0)
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        glBindRenderbuffer(GL_RENDERBUFFER, 0)

    def clearFBO(self):
        """clears old FBO"""
        if glIsFramebuffer(self.fboHandle):
            glDeleteFramebuffers(int(self.fboHandle))
        if glIsTexture(self.texHandle):
            glDeleteTextures(int(self.texHandle))

    def close(self):
        """call this to free up OpenGL resources"""
        glBindTexture(GL_TEXTURE_2D, 0)
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        glBindRenderbuffer(GL_RENDERBUFFER, 0)
        if glIsFramebuffer(self.fboHandle):
            glDeleteFramebuffers(int(self.fboHandle))
        if glIsTexture(self.texHandle):
            glDeleteTextures(int(self.texHandle))


# ===== raycast =====
rayVS = """
#version 330 core

layout(location = 1) in vec3 cubePos;
layout(location = 2) in vec3 cubeCol;

uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;

out vec4 vColor;

void main()
{
    // set position
    gl_Position = uPMatrix * uMVMatrix * vec4(cubePos.xyz, 1.0);
    // set color
    vColor = vec4(cubeCol.rgb, 1.0);
}
"""
rayFS = """
#version 330 core

in vec4 vColor;

uniform sampler2D texBackFaces;
uniform sampler3D texVolume;
uniform vec2 uWinDims;

out vec4 fragColor;

void main()
{
    // start of ray
    vec3 start = vColor.rgb;

    // calculate texture coords at fragment, which is a
    // fraction of window coords
    vec2 texc = gl_FragCoord.xy/uWinDims.xy;

    // get end of ray = back-face color
    vec3 end = texture(texBackFaces, texc).rgb;

    // calculate ray direction
    vec3 dir = end - start;

    // normalized ray direction
    vec3 norm_dir = normalize(dir);

    // the length from front to back is calculated and
    // used to terminate the ray
    float len = length(dir.xyz);

    // ray step size
    float stepSize = 0.01;

    // X-Ray projection
    vec4 dst = vec4(0.0);

    // step through the ray
    for(float t = 0.0; t < len; t += stepSize) {
        // set position to end point of ray
        vec3 samplePos = start + t*norm_dir;
        // get texture value at position
        float val = texture(texVolume, samplePos).r;
        vec4 src = vec4(val);
        // set opacity
        src.a *= 0.1;
        src.rgb *= src.a;
        // blend with previous value
        dst = (1.0 - dst.a)*src + dst;
        // exit loop when alpha exceeds threshold
        if(dst.a >= 0.95)
            break;
    }

    // set fragment color
    fragColor = dst;
}
"""


class Camera:
    """helper class for viewing"""
    def __init__(self):
        self.r = 1.5
        self.theta = 0
        self.center = [0.5, 0.5, 0.5]
        self.eye = [0.5 + self.r, 0.5, 0.5]
        self.up = [0.0, 0.0, 1.0]

    def rotate(self, clockWise):
        """rotate eye by one step"""
        if clockWise:
            self.theta = (self.theta + 5) % 360
        else:
            self.theta = (self.theta - 5) % 360
        self.eye = [0.5 + self.r*math.cos(math.radians(self.theta)),
                    0.5 + self.r*math.sin(math.radians(self.theta)),
                    0.5]


class RayCastRender:
    """class that does Ray Casting"""

    def __init__(self, width, height, volume):
        self.raycube = RayCube(width, height)
        self.width = width
        self.height = height
        self.aspect = width/float(height)
        self.program = loadShaders(rayVS, rayFS)
        self.texVolume, self.Nx, self.Ny, self.Nz = volume
        self.camera = Camera()

    def draw(self):
        pMatrix = perspective(45.0, self.aspect, 0.1, 100.0)
        mvMatrix = lookAt(self.camera.eye, self.camera.center, self.camera.up)

        # 第一遍：背面 → FBO 纹理
        texture = self.raycube.renderBackFace(pMatrix, mvMatrix)

        glUseProgram(self.program)
        glUniform2f(glGetUniformLocation(self.program, b"uWinDims"),
                    float(self.width), float(self.height))

        # 纹理单元 0：背面；单元 1：体数据
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, texture)
        glUniform1i(glGetUniformLocation(self.program, b"texBackFaces"), 0)

        glActiveTexture(GL_TEXTURE1)
        glBindTexture(GL_TEXTURE_3D, self.texVolume)
        glUniform1i(glGetUniformLocation(self.program, b"texVolume"), 1)

        # 第二遍：正面画到屏幕，片段着色器在此逐像素投射
        self.raycube.renderFrontFace(pMatrix, mvMatrix, self.program)

    def keyPressed(self, key):
        if key == 'l':
            self.camera.rotate(True)
        elif key == 'r':
            self.camera.rotate(False)

    def reshape(self, width, height):
        self.width = width
        self.height = height
        self.aspect = width/float(height)
        self.raycube.reshape(width, height)

    def close(self):
        self.raycube.close()


# ===== slicerender =====
sliceVS = """
#version 330 core

in vec3 aVert;

uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;

uniform float uSliceFrac;
uniform int uSliceMode;

out vec3 texcoord;

void main() {
  // X-slice?
  if (uSliceMode == 0) {
    texcoord = vec3(uSliceFrac, aVert.x, 1.0-aVert.y);
  }
  // Y-slice?
  else if (uSliceMode == 1) {
    texcoord = vec3(aVert.x, uSliceFrac, 1.0-aVert.y);
  }
  // Z-slice
  else {
    texcoord = vec3(aVert.x, 1.0-aVert.y, uSliceFrac);
  }

  // calculate transformed vertex
  gl_Position = uPMatrix * uMVMatrix * vec4(aVert, 1.0);
}
"""
sliceFS = """
#version 330 core

in vec3 texcoord;

uniform sampler3D tex;

out vec4 fragColor;

void main() {
  // look up color in texture
  vec4 col = texture(tex, texcoord);
  fragColor = col.rrra;
}
"""


class SliceRender:
    # slice modes
    XSLICE, YSLICE, ZSLICE = 0, 1, 2

    def __init__(self, width, height, volume):
        self.width = width
        self.height = height
        self.aspect = width/float(height)

        self.mode = SliceRender.ZSLICE

        self.program = loadShaders(sliceVS, sliceFS)
        glUseProgram(self.program)

        self.pMatrixUniform = glGetUniformLocation(self.program, b'uPMatrix')
        self.mvMatrixUniform = glGetUniformLocation(self.program, b"uMVMatrix")
        self.vertIndex = glGetAttribLocation(self.program, b"aVert")

        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)

        vertexData = numpy.array([0.0, 1.0, 0.0,
                                  0.0, 0.0, 0.0,
                                  1.0, 1.0, 0.0,
                                  1.0, 0.0, 0.0], dtype=numpy.float32)
        self.vertexBuffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glBufferData(GL_ARRAY_BUFFER, 4*len(vertexData), vertexData, GL_STATIC_DRAW)
        glEnableVertexAttribArray(self.vertIndex)
        glBindBuffer(GL_ARRAY_BUFFER, self.vertexBuffer)
        glVertexAttribPointer(self.vertIndex, 3, GL_FLOAT, GL_FALSE, 0, None)
        glBindVertexArray(0)

        self.texture, self.Nx, self.Ny, self.Nz = volume
        self.currSliceIndex = int(self.Nz/2)
        self.currSliceMax = self.Nz

    def reshape(self, width, height):
        self.width = width
        self.height = height
        self.aspect = width/float(height)

    def draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        pMatrix = ortho(-0.6, 0.6, -0.6, 0.6, 0.1, 100.0)
        # 平移 (-0.5, -0.5, -1.0)：行写平移、GL_FALSE 传入恰为标准平移矩阵
        mvMatrix = numpy.array([1.0, 0.0, 0.0, 0.0,
                                0.0, 1.0, 0.0, 0.0,
                                0.0, 0.0, 1.0, 0.0,
                                -0.5, -0.5, -1.0, 1.0], dtype=numpy.float32)
        glUseProgram(self.program)
        glUniformMatrix4fv(self.pMatrixUniform, 1, GL_FALSE, pMatrix)
        glUniformMatrix4fv(self.mvMatrixUniform, 1, GL_FALSE, mvMatrix)
        glUniform1f(glGetUniformLocation(self.program, b"uSliceFrac"),
                    float(self.currSliceIndex)/float(self.currSliceMax))
        glUniform1i(glGetUniformLocation(self.program, b"uSliceMode"),
                    self.mode)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_3D, self.texture)
        glUniform1i(glGetUniformLocation(self.program, b"tex"), 0)
        glBindVertexArray(self.vao)
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4)
        glBindVertexArray(0)

    def keyPressed(self, key):
        """key press handler"""
        if key == 'x':
            self.mode = SliceRender.XSLICE
            self.currSliceIndex = int(self.Nx/2)
            self.currSliceMax = self.Nx
        elif key == 'y':
            self.mode = SliceRender.YSLICE
            self.currSliceIndex = int(self.Ny/2)
            self.currSliceMax = self.Ny
        elif key == 'z':
            self.mode = SliceRender.ZSLICE
            self.currSliceIndex = int(self.Nz/2)
            self.currSliceMax = self.Nz
        elif key == 'l':  # 左箭头：索引 +1（原书正文说反了，以代码为准）
            self.currSliceIndex = (self.currSliceIndex + 1) % self.currSliceMax
        elif key == 'r':
            self.currSliceIndex = (self.currSliceIndex - 1) % self.currSliceMax

    def close(self):
        pass


# ===== volrender =====
class RenderWin:
    """GLFW Rendering window class"""
    def __init__(self, imageDir):
        cwd = os.getcwd()  # glfwInit 会改变当前工作目录，先保存
        glfw.glfwInit()
        os.chdir(cwd)

        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MAJOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_CONTEXT_VERSION_MINOR, 3)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE)
        glfw.glfwWindowHint(glfw.GLFW_OPENGL_PROFILE, glfw.GLFW_OPENGL_CORE_PROFILE)

        self.width, self.height = 512, 512
        self.aspect = self.width/float(self.height)
        self.win = glfw.glfwCreateWindow(self.width, self.height, b"volrender")
        glfw.glfwMakeContextCurrent(self.win)

        glViewport(0, 0, self.width, self.height)
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.0, 0.0, 0.0, 0.0)

        glfw.glfwSetMouseButtonCallback(self.win, self.onMouseButton)
        glfw.glfwSetKeyCallback(self.win, self.onKeyboard)
        glfw.glfwSetWindowSizeCallback(self.win, self.onSize)

        self.volume = loadVolume(imageDir)
        self.renderer = RayCastRender(self.width, self.height, self.volume)
        self.exitNow = False

    def onMouseButton(self, win, button, action, mods):
        pass

    def onKeyboard(self, win, key, scancode, action, mods):
        if key == glfw.GLFW_KEY_ESCAPE:  # 原书用 is 比较常量，应为 ==
            self.renderer.close()
            self.exitNow = True
        else:
            if action == glfw.GLFW_PRESS or action == glfw.GLFW_REPEAT:
                if key == glfw.GLFW_KEY_V:
                    # toggle render mode
                    if isinstance(self.renderer, RayCastRender):
                        self.renderer = SliceRender(self.width, self.height,
                                                    self.volume)
                    else:
                        self.renderer = RayCastRender(self.width, self.height,
                                                      self.volume)
                    self.renderer.reshape(self.width, self.height)
                else:
                    # 键码转字符：GLFW 依赖只留在主文件
                    keyDict = {glfw.GLFW_KEY_X: 'x', glfw.GLFW_KEY_Y: 'y',
                               glfw.GLFW_KEY_Z: 'z',
                               glfw.GLFW_KEY_LEFT: 'l', glfw.GLFW_KEY_RIGHT: 'r'}
                    try:
                        self.renderer.keyPressed(keyDict[key])
                    except KeyError:
                        pass

    def onSize(self, win, width, height):
        self.width, self.height = width, height
        self.aspect = width/float(height)
        glViewport(0, 0, self.width, self.height)
        self.renderer.reshape(width, height)

    def run(self):
        while not glfw.glfwWindowShouldClose(self.win) and not self.exitNow:
            self.renderer.draw()
            glfw.glfwSwapBuffers(self.win)
            glfw.glfwWaitEvents()  # 事件驱动：无动画，等事件再重画
        glfw.glfwTerminate()


def main():
    print('starting volrender...')
    parser = argparse.ArgumentParser(description="Volume Rendering...")
    parser.add_argument('--dir', dest='imageDir', required=True)
    args = parser.parse_args()
    rwin = RenderWin(args.imageDir)
    rwin.run()


if __name__ == '__main__':
    main()
