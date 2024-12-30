////////////////////////////////////////////////////////////////////////
//
//   KAIST, Spring 2023
//   CS380: Introduction to Computer Graphics
//   Instructor: Minhyuk Sung (mhsung@kaist.ac.kr)
//   Last Update: Juil Koo (63days@kaist.ac.kr)
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
// If your OS is LINUX, uncomment the line below.
//#include <tr1/memory>

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
// #else
// #   include <GL/glut.h>
#endif

#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "quat.h"
#include "rigtform.h"
#include "arcball.h"

#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"

using namespace std;      // for string, vector, iostream, and other standard C++ stuff
// If your OS is LINUX, uncomment the line below.
//using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------
const bool g_Gl2Compatible = true;


static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

static const int PICKING_SHADER = 2; // index of the picking shader is g_shaerFiles
static const int g_numShaders = 3; // 3 shaders
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
  }
};

// a typedef that declares
// our own Shape node which draws using our Geometry.
typedef SgGeometryShapeNode<Geometry> MyShapeNode;

// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere;

// Declare the scene graph and pointers to suitable nodes in the scene graph
static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking


// --------- Scene

static const int numCubes = 2;
static const int numViews = 3;
static int idxView = 0; // indexing for view changing when 'v' is pressed
static int idxFrame = 0; // indexing for object manipulating when 'm' is pressed
static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static Cvec3f g_objectColors[3] = {Cvec3f(1, 0, 0), Cvec3f(0, 1, 0), Cvec3f(0, 1, 1)}; // colors of objs
static RigTForm g_eyeRbt;
static RigTForm g_sphRbt = RigTForm(Cvec3(0,0,0));
static double g_arcballScreenRadius = 128.0, g_arcballScale = 1.0;
static bool g_isPicking = false;
static bool g_isArcball = true;


enum CAMERA_STATUS {
  SKY_CAMERA,
  CUBE1_CAMERA,
  CUBE2_CAMERA,
};

enum WORLDSKYFRAME_STATUS {
  WORLD_SKY,
  SKY_SKY,
};

///////////////// END OF G L O B A L S //////////////////////////////////////////////////

static shared_ptr<SgRbtNode> eyeNode() {
  if (idxView == SKY_CAMERA) {
    return g_skyNode;
  }
  else if (idxView == CUBE1_CAMERA) {
    return g_robot1Node;
  }
  else { // CUBE2_CAMERA
    return g_robot2Node;
  }
}


static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSpheres() {
  int slices = 20, stacks = 20; // input
  int ibLen, vbLen; // output
  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeSphere(1, slices, stacks, vtx.begin(), idx.begin());
  g_sphere.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS380_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff(const ShaderState& curSS, bool picking) {
  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // eyeRbt
  RigTForm eyeRbt = getPathAccumRbt(g_world, eyeNode());
  g_eyeRbt = eyeRbt;

  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  // Replace the code for drawing the ground the two cubes with:
  if (!picking) {
    Drawer drawer(invEyeRbt, curSS);
    g_world->accept(drawer);

    const RigTForm groundRbt = RigTForm();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // draw wireframe

    // if manipulating a cube but is not the current eyeframe
    if (g_currentPickedRbtNode != shared_ptr<SgRbtNode>()) {
      RigTForm MVM = invEyeRbt * getPathAccumRbt(g_world, g_currentPickedRbtNode);
      Matrix4 scale = Matrix4::makeScale(g_arcballScale * g_arcballScreenRadius);
      Matrix4 NMVM = normalMatrix(rigTFormToMatrix(MVM) * scale); // radius = g_arcballScale * g_arcballScreenRadius;
      sendModelViewNormalMatrix(curSS, rigTFormToMatrix(MVM) * scale, NMVM);
      safe_glUniform3f(curSS.h_uColor, g_objectColors[2][0], g_objectColors[2][1], g_objectColors[2][2]);
      g_sphere->draw(curSS);
      g_isArcball = true;
    }

    else if (g_isArcball) {
      RigTForm MVM = invEyeRbt * groundRbt;
      if (!g_mouseMClickButton && !(g_mouseLClickButton && g_mouseRClickButton))
        g_arcballScale = getScreenToEyeScale(MVM.getTranslation()[2], g_frustFovY, g_windowHeight);
      Matrix4 scale = Matrix4::makeScale(g_arcballScale * g_arcballScreenRadius);
      Matrix4 NMVM = normalMatrix(rigTFormToMatrix(MVM) * scale); // radius = g_arcballScale * g_arcballScreenRadius;
      sendModelViewNormalMatrix(curSS, rigTFormToMatrix(MVM) * scale, NMVM);
      safe_glUniform3f(curSS.h_uColor, g_objectColors[2][0], g_objectColors[2][1], g_objectColors[2][2]);
      g_sphere->draw(curSS);
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // draw filled again
  }
  else {
    cout << "picked!" << endl;
    Picker picker(invEyeRbt, curSS);
    g_world->accept(picker);
    glFlush();
    g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
    if (g_currentPickedRbtNode == g_groundNode) {
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>();   // set to NULL
      g_isArcball = false;
    }
  }
  
}


static void pick() {
  // We need to set the clear color to black, for pick rendering.
  // so let's save the clear color
  GLdouble clearColor[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);

  glClearColor(0, 0, 0, 0);

  // using PICKING_SHADER as the shader
  glUseProgram(g_shaderStates[PICKING_SHADER]->program);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff(*g_shaderStates[PICKING_SHADER], true);

  // Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
  // to see result of the pick rendering pass
  // glutSwapBuffers();

  //Now set back the clear color
  glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);

  checkGlErrors();
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(*g_shaderStates[g_activeShader], false);

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;

  g_arcballScreenRadius = 0.25 * min(g_windowHeight, g_windowWidth);

  updateFrustFovY();
  glutPostRedisplay();
}

Cvec3 isIntersect(double x, double y) {
  double t1, t2, t;
  double radius = g_arcballScale * g_arcballScreenRadius;
  
  Cvec3 d = getModelViewRay(Cvec2(x, y), g_frustFovY, g_windowWidth, g_windowHeight);
  Cvec3 sph_origin = g_sphRbt.getTranslation();
  Cvec3 ray_origin = g_eyeRbt.getTranslation();
  Cvec3 oc = sph_origin - ray_origin;
  oc = Cvec3(inv(g_eyeRbt).getRotation() * Cvec4(oc, 0));
  double a = dot(d, d), h = dot(d, oc);
  double c = dot(oc, oc) - radius * radius;
  double discriminant = h*h - a*c;

  if (discriminant >= 0) {
    t1 = (h + sqrt(discriminant)) / a;
    t2 = (h - sqrt(discriminant)) / a;

    t = fabs(t1) < fabs(t2) ? t1 : t2;
    Cvec3 v = - oc + d * t; // ray_origin + d * t - sph_origin
    return v;
  }
  else {
    return normalize(d * dot(oc, d) - oc) * radius; // return closest point
  }
}

static void motion(const int x, const int y) {
  const double dx = x - g_mouseClickX, x1 = g_mouseClickX, x2 = x;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY, y1 = g_mouseClickY, y2 = g_windowHeight - y - 1;

  g_eyeRbt = getPathAccumRbt(g_world, eyeNode()); // must to do: g_isarcball 설정 -> picking, world-sky 시 arcball draw + rot/trans되도록 설정

  if (g_isArcball && (g_currentPickedRbtNode != shared_ptr<SgRbtNode>())) {
    g_sphRbt = getPathAccumRbt(g_world, g_currentPickedRbtNode);
  }
  else if (idxView == SKY_CAMERA && idxFrame == SKY_SKY) {
    g_sphRbt = g_eyeRbt;
  }
  else {
    g_sphRbt = RigTForm();
  }

  RigTForm m, eyeRbt, sphRbt, tmp;
  if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
    if (g_isArcball) { // calculate v1, v2 from projection infos. 
      // first method: Quadratic solution
      Cvec3 v1 = isIntersect(x1, y1), v2 = isIntersect(x2, y2);

      // second method: Geometric solution
      // RigTForm MVM = inv(g_eyeRbt) * g_sphRbt;
      // Cvec2 sphPts = getScreenSpaceCoord(MVM.getTranslation(), makeProjectionMatrix(), g_frustNear, g_frustFovY, g_windowWidth, g_windowHeight);
      // double z1 = sqrt(max(0.0, g_arcballScreenRadius * g_arcballScreenRadius - (x1 - sphPts[0]) * (x1 - sphPts[0]) - (y1 - sphPts[1]) * (y1 - sphPts[1])));
      // double z2 = sqrt(max(0.0, g_arcballScreenRadius * g_arcballScreenRadius - (x2 - sphPts[0]) * (x2 - sphPts[0]) - (y2 - sphPts[1]) * (y2 - sphPts[1])));
      // Cvec3 v1 = normalize(Cvec3(x1 - sphPts[0], y1 - sphPts[1], z1));
      // Cvec3 v2 = normalize(Cvec3(x2 - sphPts[0], y2 - sphPts[1], z2));

      m = RigTForm(Quat(0, normalize(v2)) * Quat(0, normalize(v1 * -1)));
    }
    else
      m = RigTForm(Quat::makeXRotation(-dy) * Quat::makeYRotation(dx));
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if (g_isArcball) m = RigTForm(Cvec3(dx, dy, 0) * g_arcballScale);
    else m = RigTForm(Cvec3(dx, dy, 0) * 0.01);
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    if (g_isArcball) m = RigTForm(Cvec3(0, 0, -dy) * g_arcballScale);
    else m = RigTForm(Cvec3(0, 0, -dy) * 0.01);
  }

  if (g_mouseClickDown) {
    RigTForm a;
    if (g_isArcball) {
      if (idxView == SKY_CAMERA && idxFrame == WORLD_SKY && (g_currentPickedRbtNode == shared_ptr<SgRbtNode>())) { 
          RigTForm orgRbt;
          a = makeMixedFrame(orgRbt, g_eyeRbt);
          tmp = doMtoOwrtA(inv(m), g_eyeRbt, a);
          eyeNode()->setRbt(tmp);
        }
      else { // move the object
        a = makeMixedFrame(getPathAccumRbt(g_world, g_currentPickedRbtNode), g_eyeRbt);
        tmp = doMtoOwrtA(m, g_currentPickedRbtNode->getRbt(), inv(getPathAccumRbt(g_world, g_currentPickedRbtNode, 1)) * a);
        g_currentPickedRbtNode->setRbt(tmp);
      }
    }
    else {
      // move the eye
      tmp = doMtoOwrtA(m, g_eyeRbt, g_eyeRbt);
      eyeNode()->setRbt(tmp);
    }
    
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;

  if (g_isPicking) {
    pick();
    g_isPicking = false;
  }
  glutPostRedisplay(); // we always redraw if we changed the scene
}


static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "v\t\tCycle view\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case 'v':
    idxView = (idxView + 1) % 3;
    if (idxView == SKY_CAMERA) {
      cout << "Sky Camera" << endl;
    }
    else if (idxView == CUBE1_CAMERA) {
      cout << "Robot1 Camera" << endl;
    }
    else { // CUBE2_CAMERA
      cout << "Robot2 Camera" << endl;
    }
    break;
  case 'm':
    if (idxView == SKY_CAMERA) {
      idxFrame = (idxFrame + 1) % 2;
      if (idxFrame == SKY_SKY) {
        g_isArcball = false;
        cout << "Sky-sky frame" << endl;
      } else { // WORLD_SKY
        g_isArcball = true;
        cout << "World-sky frame" << endl;
      }
    } else {
      cout << "'m' must be used when the view and the object are all sky-camera" << endl;
    }
    break;
  case 'p':
    g_isPicking = true;
    cout << "Picking" << endl;
    break;
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("CS380: Assignment 2");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes();
  initSpheres();
}

static void constructRobot(shared_ptr<SgTransformNode> base, const Cvec3& color) {
  const double ARM_LEN = 0.7,
               ARM_THICK = 0.25,
               TORSO_LEN = 1.5,
               TORSO_THICK = 0.25,
               TORSO_WIDTH = 1,
               HEAD_RADIUS = 0.5;
  const int NUM_JOINTS = 10,
            NUM_SHAPES = 10;

  struct JointDesc {
    int parent;
    float x, y, z;
  };

  JointDesc jointDesc[NUM_JOINTS] = {
    {-1}, // torso
    {0,  TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper right arm
    {1,  ARM_LEN, 0, 0}, // lower right arm
    {0,  -TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper right arm
    {3,  -ARM_LEN, 0, 0}, // lower right arm
    {0, (TORSO_WIDTH) / 4, -(TORSO_LEN / 2), 0}, // upper right leg
    {5, 0, -ARM_LEN, 0},  // lower right leg
    {0, -(TORSO_WIDTH) / 4, -(TORSO_LEN / 2), 0},    // upper left leg
    {7, 0, -ARM_LEN, 0},        // lower left leg
    {0, 0, (TORSO_LEN / 2), 0}  // head
  };

  struct ShapeDesc {
    int parentJointId;
    float x, y, z, sx, sy, sz;
    shared_ptr<Geometry> geometry;
  };

  ShapeDesc shapeDesc[NUM_SHAPES] = {
    {0, 0,         0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube}, // torso
    {1, ARM_LEN/2, 0, 0, ARM_LEN / 2, ARM_THICK / 2, ARM_THICK / 2, g_sphere}, // upper right arm
    {2, ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // lower right arm
    {3, -ARM_LEN/2, 0, 0, ARM_LEN / 2, ARM_THICK / 2, ARM_THICK / 2, g_sphere}, // upper left arm
    {4, -ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // lower left arm
    {5, 0, -ARM_LEN/2, 0, ARM_THICK / 2, ARM_LEN / 2, ARM_THICK / 2, g_sphere}, // upper right leg
    {6, 0, -ARM_LEN/2, 0, ARM_THICK, ARM_LEN, ARM_THICK, g_cube}, // lower right leg
    {7, 0, -ARM_LEN/2, 0, ARM_THICK / 2, ARM_LEN / 2, ARM_THICK / 2, g_sphere}, // upper left leg
    {8, 0, -ARM_LEN/2, 0, ARM_THICK, ARM_LEN, ARM_THICK, g_cube}, // lower left leg
    {9, 0, HEAD_RADIUS, 0, HEAD_RADIUS, HEAD_RADIUS, HEAD_RADIUS, g_sphere},
  };

  shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];

  for (int i = 0; i < NUM_JOINTS; ++i) {
    if (jointDesc[i].parent == -1)
      jointNodes[i] = base;
    else {
      jointNodes[i].reset(new SgRbtNode(RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
      jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
    }
  }
  for (int i = 0; i < NUM_SHAPES; ++i) {
    shared_ptr<MyShapeNode> shape(
      new MyShapeNode(shapeDesc[i].geometry,
                      color,
                      Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z),
                      Cvec3(0, 0, 0),
                      Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
    jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
  }
}

static void initScene() {
  g_world.reset(new SgRootNode());

  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 6.0))));

  g_groundNode.reset(new SgRbtNode());
  g_groundNode->addChild(shared_ptr<MyShapeNode>(
                           new MyShapeNode(g_ground, Cvec3(0.1, 0.95, 0.1))));

  g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
  g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));

  constructRobot(g_robot1Node, Cvec3(1, 0, 0)); // a Red robot
  constructRobot(g_robot2Node, Cvec3(0, 0, 1)); // a Blue robot

  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_robot1Node);
  g_world->addChild(g_robot2Node);
}


int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    glewInit(); // load the OpenGL extensions

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

    initGLState();
    initShaders();
    initGeometry();
    initScene();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
