// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>

// Linear Algebra Library
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h> // use M_PI

#include <iostream>

#include <vector>

// Timer
#include <chrono>

using namespace std;
using namespace Eigen;

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_C;

// Vertex positions initialization
Eigen::MatrixXf V(2,6);

// Contains the view transformation
Eigen::Matrix4f view(4,4);

// MatrixXf Models(4,4);
MatrixXf identityMatrix(4,4);

// Contains the per-vertex color
Eigen::MatrixXf C(3,3);
Eigen::Matrix3f originalColor(3,3); 

// Initiate total number of triangle
int num = 0;

// Insertion
bool insertionActive = false;
int insertionCount = 0; // The nth click in triangle insertion

// Translation
bool translationActive = false;
bool translationClick = false;
bool translationSelected = false;
double currentX, currentY, previousX, previousY;

// Deletion
bool deletionActive = false;

// Rotation/scale
bool clockwise = false;
bool counterclockwise = false;
bool scaleUp = false;
bool scaleDown = false;

bool colorActive = false;
int c_vertex = -1; // Color the nth vertice
// Number selected in the color mode
int selectedColor = -1;
MatrixXf colors(3, 9);
// Initialize selected vertices
int vertex_click1 = -1;
int vertex_click2 = -1;
int vertex_click3 = -1;
int vertex_delete1 = -1;
int vertex_delete2 = -1;
int vertex_delete3 = -1;

// starting point of the model matrix for the 1st vertex 
int vertexModel = -1;

// Animation
bool animationActive = false;
bool animationClicked = false;
bool animate = false;
double beginX = -1; double beginY = -1;
double endX = -1; double endY = -1;

// initialize starting time 
auto t_start = std::chrono::high_resolution_clock::now();

// Translation matrices initialization
Eigen::MatrixXf translation(4,4);
Eigen::MatrixXf rotation(4,4);
Eigen::MatrixXf scaling(4,4);
Eigen::MatrixXf model(4,4);

// Previous translation matrix for an animated triangle
Eigen::Matrix4f original(4,4);

// Map an index in V to the corresponding index in model
int VtoModel(int v_index){
  return v_index + floor(v_index/3);
}

// Linear interpolation between keyframes
float linearInterpolate(float a1,float a2, float mu){
    return (a1 * (1-mu) + a2 * mu);
}

void resetTranslationVariables() // Reset the vertices number in translation mode
{
  if(vertex_click1 > -1)
  {
    C.col(vertex_click1) = originalColor.col(0);
    C.col(vertex_click2) = originalColor.col(1);
    C.col(vertex_click3) = originalColor.col(2);
    if(!animationActive){
      vertex_click1 = -1;
      vertex_click2 = -1;
      vertex_click3 = -1;
      vertexModel = -1;
    }
    translationSelected = false;
    VBO_C.update(C);
  }
}

void translateFun()
{
  // translation based on previous and current coordinates
  float x_diff = currentX - previousX;
  float y_diff = currentY - previousY;
  translation(0,vertexModel+3) = translation(0,vertexModel+3) + x_diff;
  translation(1,vertexModel+3) = translation(1,vertexModel+3) + y_diff;
  // Eigen block: Block of size (p,q), starting at (i,j) ----  matrix.block(i,j,p,q);
  model.block(0, vertexModel, 4,4) = translation.block(0, vertexModel, 4, 4) * rotation.block(0, vertexModel, 4, 4) * scaling.block(0, vertexModel, 4, 4);
}

void animateTriangle()
{
  auto t_now = std::chrono::high_resolution_clock::now();
  float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
  if(animationActive){
    if(time <= 1.0) //Restrict the range
    {
      float x_point =  linearInterpolate(beginX, endX, time);
      float y_point =  linearInterpolate(beginY, endY, time);

      // Change V pos to model position and do the translation
      Vector4f V_pos(V(0, vertex_click1), V(1, vertex_click1), 0.0, 1.0);
      Vector4f pre_pos = model.block(0, vertexModel, 4, 4) * V_pos;

      float x_diff = x_point - pre_pos[0];
      float y_diff = y_point - pre_pos[1];
      translation(0,vertexModel+3) += x_diff;
      translation(1,vertexModel+3) += y_diff;
      // Update model
      model.block(0, vertexModel, 4,4) = translation.block(0, vertexModel, 4, 4) * rotation.block(0, vertexModel, 4, 4) * scaling.block(0, vertexModel, 4, 4);}
      else
      {
      animationActive = false;
      animationClicked = false;
      animate = false;
      resetTranslationVariables();
    }
  }
}

void panFun(int dir,int windowWidth, int windowHeight)
{
  //Vertical pan
  Vector4f p_screen(0.0, (0.2 * windowHeight), 0.0, 0.0);
  Vector4f p_canonical( ((p_screen[0]/windowWidth) * 2 -1), (p_screen[1]/windowHeight)*2-1, 0.0, 0.0);
  Vector4f p_world = view.inverse() * p_canonical;
  //Horizontal pan
  Vector4f p_screen2((0.2 * windowWidth), 0.0, 0.0, 0.0);
  Vector4f p_canonical2( ((p_screen2[0]/windowWidth) * 2 -1), (p_screen2[1]/windowHeight)*2-1, 0.0, 0.0);
  Vector4f p_world2 = view.inverse() * p_canonical2;
  if(dir == 1){
    view(1,3) = view(1,3) - p_world[1];}
  else if(dir == 2){
    view(0,3) = view(0,3) + p_world2[0];}
  else if(dir == 3){
    view(1,3) = view(1,3) + p_world[1];}
  else if(dir == 4){
    view(0,3) =view(0,3) - p_world2[0];}
}

void zoom(int z, int windowWidth, int windowHeight)
{
  if(z == 1)
  {
    for(int i = 0; i<3; i++)
    {view(i,i) *= 1.2;}
  }
  else
  {
    for(int i = 0; i<3; i++)
    {view(i,i) *= 0.8;}
  }
}

void colorVertex()
{
  if(colorActive && c_vertex > -1)
  {
    int idx = selectedColor - 1 ; //offset by 1 to reach the entries in col 0
    C.col(c_vertex) = colors.col(idx);
    VBO_C.update(C);
  }
}

// Calculate the center coordinates of the triangle
Vector2f barycenterCalc()
{
  // Calculate center
  float barycenterX = (V(0, vertex_click1)+V(0, vertex_click2)+V(0, vertex_click3)) / 3;
  float barycenterY = (V(1, vertex_click1)+V(1, vertex_click2)+V(1, vertex_click3)) / 3;
  Vector2f center(barycenterX, barycenterY);
  return center;
}

//https://open.gl/transformations
void scaleFun(bool scaleUp)
{
  if(translationActive && translationSelected)
  {
    Vector2f c = barycenterCalc();
    float s;
    if(scaleUp)
    {s = 1.25;}
    else
    {s = 0.75;}
    Matrix4f toScale;
    toScale <<
    s, 0, 0, 0,
    0, s, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    // Translate the barycenter to ensure the scaling happens around the center of the triangle
    Matrix4f movetoO;
    movetoO <<
    1, 0, 0, -c[0],
    0, 1, 0, -c[1],
    0, 0, 1, 0,
    0, 0, 0, 1;

    Matrix4f moveBack;
    moveBack <<
    1, 0, 0, c[0],
    0, 1, 0, c[1],
    0, 0, 1, 0,
    0, 0, 0, 1;

    scaling.block(0, vertexModel, 4, 4) = moveBack * toScale * movetoO * scaling.block(0, vertexModel, 4, 4);
    model.block(0, vertexModel, 4,4) = translation.block(0, vertexModel, 4, 4) * rotation.block(0, vertexModel, 4, 4) * scaling.block(0, vertexModel, 4, 4) ;
  }
}

void rotateFun(bool clockwise)
{
  if(translationActive && translationSelected)
  {
    Vector2f c = barycenterCalc();
    float theta;
    if(clockwise)
    {theta = -10 * (M_PI / 180);}
    else
    {theta = 10 * (M_PI / 180);}
    Matrix4f rotate;
    rotate <<
    cos(theta), -sin(theta),  0, 0,
    sin(theta), cos(theta),   0, 0,
    0,          0,            1, 0,
    0,          0,            0, 1;

    Matrix4f movetoO; // move to the origin
    movetoO <<
    1, 0, 0, -c[0],
    0, 1, 0, -c[1],
    0, 0, 1, 0,
    0, 0, 0, 1;

    Matrix4f moveBack; // move back
    moveBack <<
    1, 0, 0, c[0],
    0, 1, 0, c[1],
    0, 0, 1, 0,
    0, 0, 0, 1;

    rotation.block(0, vertexModel, 4, 4) = moveBack * rotate * movetoO * rotation.block(0, vertexModel, 4, 4);
    model.block(0, vertexModel, 4,4) = translation.block(0, vertexModel, 4, 4) * rotation.block(0, vertexModel, 4, 4) * scaling.block(0, vertexModel, 4, 4) ;

  }
}

// Judge if the mouse click happens on a triangle
bool clickOnTriangle(double mousex, double mousey, float coord1_x, float coord1_y, float coord2_x, float coord2_y, float coord3_x, float coord3_y)
{
  Matrix3f A;
  Vector3f b;
  A << coord1_x, coord2_x, coord3_x, coord1_y, coord2_y, coord3_y, 1, 1, 1;
  b << mousex, mousey, 1;
  Vector3f x = A.colPivHouseholderQr().solve(b);
  float alpha = x[0];
  float beta = x[1];
  float gamma = x[2];
  if(alpha > 0 && beta > 0 && gamma > 0)
  {return true;}
  else
  {return false;}
}

void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
  // Get the size of the window
  int width, height;
  glfwGetWindowSize(window, &width, &height);

  // Convert screen position to world coordinates
  Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
  Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
  Eigen::Vector4f p_world = view.inverse()*p_canonical;

  double xworld = p_world[0];
  double yworld = p_world[1];

  // Track the mouse positions
  if(!previousX && !previousY)
  {
    previousX = xworld;
    previousY = yworld;
  }
  else
  {
    previousX = currentX;
    previousY = currentY;
    currentX = xworld;
    currentY = yworld;
  }
  if(insertionActive && insertionCount > 0)
  {
    // Store coordinates in V and send to GPU
    if(insertionCount == 1){
      V.col((num*3) + 1) << xworld, yworld;
    }else if(insertionCount == 2){ // the end of 2 sides of the triangle
      V.col((num*3) + 3) << xworld, yworld;
      V.col((num*3) + 5) << xworld, yworld;
    }
    VBO.update(V);
  }
  if(translationClick||animationClicked)
  translateFun();
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // Convert screen position to world coordinates
    Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
    Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    Eigen::Vector4f p_world = view.inverse()*p_canonical;

    double xworld = p_world[0];
    double yworld = p_world[1];

    if(insertionActive && action == GLFW_RELEASE) // creating a triangle by 6 points, 3 lines
    {
      if(insertionCount == 0) // First click
      {
        V.conservativeResize(2, (num*3) + 6); // use conservative resize to make sure the original entries do not change
        C.conservativeResize(3, (num*3) + 6);

        // Show white lines for insertion tracking
        for (int i = 0; i < 6; i++)
        {
          C.col(num*3+i) << 1.0, 1.0, 1.0;
        }

        // 1st click, set all points to the cursor coordinate
        for (int i = 0; i < 6; i++)
        {
          V.col(num*3+i) << xworld, yworld;
        }
        VBO.update(V);
        VBO_C.update(C);
      }
      // 2nd click
      else if(insertionCount == 1){
        // Change V for multiple dynamic lines" << std::endl;
        V.col(num*3 + 1) << xworld, yworld;
        V.col(num*3 + 2) << V(0, num*3),V(1,num*3);
        V.col(num*3 + 4) << xworld, yworld;// Show the third line
        VBO.update(V);
      }
      // 3rd click
      else if(insertionCount == 2){
        //  Construct a new triangle, which requires additional space for V and C.
        Eigen::MatrixXf temp_V(2, (num*3) + 3);
        Eigen::MatrixXf temp_C(3, (num*3) + 3);

        int start = num*3;
        int model_start = num*4;

        // copy the previous V and C to the temporary matrix
        for(int i = 0; i < start; i++)
        {
          temp_V.col(i) = V.col(i);
          temp_C.col(i) = C.col(i);
        }

        temp_V.col(start) << V(0, start), V(1, start);
        temp_V.col(start + 1) << V(0, (start + 1)), V(1, (start + 1));
        temp_V.col(start + 2) << xworld, yworld;

        temp_C.col(start) << 1.0, 0.0, 0.0;
        temp_C.col(start + 1) << 1.0, 0.0, 0.0;
        temp_C.col(start + 2) << 1.0, 0.0, 0.0;

        num++;

        model.conservativeResize(4, (num*4));
        scaling.conservativeResize(4, (num*4));
        translation.conservativeResize(4, (num*4));
        rotation.conservativeResize(4, (num*4));

        model.col(model_start) << 1, 0, 0, 0;
        model.col(model_start+1) << 0, 1, 0, 0;
        model.col(model_start+2) << 0, 0, 1, 0;
        model.col(model_start+3) << 0, 0, 0, 1;

        scaling.col(model_start) << 1, 0, 0, 0;
        scaling.col(model_start+1) << 0, 1, 0, 0;
        scaling.col(model_start+2) << 0, 0, 1, 0;
        scaling.col(model_start+3) << 0, 0, 0, 1;

        translation.col(model_start) << 1, 0, 0, 0;
        translation.col(model_start+1) << 0, 1, 0, 0;
        translation.col(model_start+2) << 0, 0, 1, 0;
        translation.col(model_start+3) << 0, 0, 0, 1;

        rotation.col(model_start) << 1, 0, 0, 0;
        rotation.col(model_start+1) << 0, 1, 0, 0;
        rotation.col(model_start+2) << 0, 0, 1, 0;
        rotation.col(model_start+3) << 0, 0, 0, 1;

        // Update V and C
        V.resize(2, (num*3));
        C.resize(3, (num*3));
        V = temp_V;
        C = temp_C;
        VBO.update(V);
        VBO_C.update(C);
      }
      insertionCount++;
    }

    else if(translationActive)
    {
      if(action == GLFW_PRESS)
      {
        // Judge if the cursor is in a triangle
        for(int i = 0; i < V.cols(); i = i+3) // For each triangle
        {
          int model_idx = VtoModel(i);

          Vector4f point1(V(0,i), V(1,i), 0, 1);
          Vector4f point2(V(0,i + 1), V(1,i + 1), 0, 1);
          Vector4f point3(V(0,i + 2), V(1,i + 2), 0, 1);

          Vector4f p1onScreen = model.block(0, model_idx, 4, 4) * point1;
          Vector4f p2onScreen = model.block(0, model_idx, 4, 4) * point2;
          Vector4f p3onScreen = model.block(0, model_idx, 4, 4) * point3;

          translationClick = clickOnTriangle(xworld, yworld, p1onScreen[0], p1onScreen[1], p2onScreen[0], p2onScreen[1], p3onScreen[0], p3onScreen[1]);

          if(translationClick)
          {
            if(vertex_click1 != i){
              if(vertex_click1 > -1){
                C.col(vertex_click1) = originalColor.col(0);
                C.col(vertex_click2) = originalColor.col(1);
                C.col(vertex_click3) = originalColor.col(2);
              }
              vertex_click1 = i;
              vertex_click2 = i+1;
              vertex_click3 = i+2;

              vertexModel = VtoModel(i);

              // Store the original color
              originalColor.col(0) = C.col(vertex_click1);
              originalColor.col(1) = C.col(vertex_click2);
              originalColor.col(2) = C.col(vertex_click3);

              // Mark the selected triangle in blue
              C.col(vertex_click1) << 0.0, 0.0, 1.0;
              C.col(vertex_click2) << 0.0, 0.0, 1.0;
              C.col(vertex_click3) << 0.0, 0.0, 1.0;
              VBO_C.update(C);
            }
            break;
          }
        }
      }
      else if(action == GLFW_RELEASE && translationClick)
      {
        translationClick = false;
        translationSelected = true;
        VBO_C.update(C);
      }
    }

    else if(deletionActive && action == GLFW_RELEASE)
    {
      // Judge if the cursor is in the triangle
      for(int i = 0; i < V.cols(); i=i+3) // 3 vertices per triangle
      {
        int model_idx = VtoModel(i);

        Vector4f point1(V(0,i), V(1,i), 0, 1);
        Vector4f point2(V(0,i + 1), V(1,i + 1), 0, 1);
        Vector4f point3(V(0,i + 2), V(1,i + 2), 0, 1);

        Vector4f p1onScreen = model.block(0, model_idx, 4, 4) * point1;
        Vector4f p2onScreen = model.block(0, model_idx, 4, 4) * point2;
        Vector4f p3onScreen = model.block(0, model_idx, 4, 4) * point3;

        // Delete if the click is in a triangle
        if (clickOnTriangle(xworld, yworld, p1onScreen[0], p1onScreen[1], p2onScreen[0], p2onScreen[1], p3onScreen[0], p3onScreen[1]))
        {
          // Remove vertices and update V
          num--;
          vertex_delete1 = i;
          vertex_delete2 = i + 1;
          vertex_delete3 = i + 2;
          
          int c = 0; // counter for temporary matrices
          int model_delete = VtoModel(i);
          Eigen::MatrixXf temp_model(4, (num*4));
          Eigen::MatrixXf temp_scaling(4, (num*4));
          Eigen::MatrixXf temp_rotation(4, (num*4));
          Eigen::MatrixXf temp_translation(4, (num*4));
          Eigen::MatrixXf temp_V(2, (num*3));
          Eigen::MatrixXf temp_C(3, (num*3));
          
          // Copy the rest entries V into temp_V
          for(int i = 0; i < V.cols(); i++)
          {
            if(i == vertex_delete1 || i == vertex_delete2 || i == vertex_delete3)
            {continue;}
            else
            {
              temp_V.col(c) = V.col(i);
              temp_C.col(c) = C.col(i);
              c++;
            }
          }
          
          c = 0;
          for(int i = 0; i < model.cols(); i++){
            if(i == model_delete || i == model_delete+1 || i == model_delete+2 || i == model_delete+3)
            {continue;}
            else
            {
              temp_model.col(c) = model.col(i);
              temp_scaling.col(c) = scaling.col(i);
              temp_rotation.col(c) = rotation.col(i);
              temp_translation.col(c) = translation.col(i);
              c++;
            }
          }
          model = temp_model;
          scaling = temp_scaling;
          rotation = temp_rotation;
          translation = temp_translation;
          V = temp_V;
          C = temp_C;
          break;
        }
      }
      VBO.update(V);
      VBO_C.update(C);
    }

    else if(colorActive && action == GLFW_RELEASE)
    {
      // Find the closest vertex
      float smallestd = 10;
      int model_idx = 0;
      for(int i = 0; i < V.cols(); i++)
      {
        if(i%3 == 0) 
        {model_idx = VtoModel(i);}
        Vector4f p(V(0,i), V(1,i), 0, 1);
        Vector4f ponScreen = model.block(0, model_idx, 4, 4) * p;
        float dX = xworld - ponScreen[0];
        float dY = yworld - ponScreen[1];
        float d = sqrt(pow(dX,2)+pow(dY,2));
        if(d < smallestd)
        {
          smallestd = d;
          c_vertex = i;
        }
      }
    }

    else if(animationActive)
    {
      if(action == GLFW_PRESS && !animationClicked){
        // Judge if the cursor is inside a triangle
        for(int i = 0; i < V.cols(); i=i+3) // 3 vertices per triangle
        {
          int model_idx = VtoModel(i);

          Vector4f point1(V(0,i), V(1,i), 0, 1);
          Vector4f point2(V(0,i + 1), V(1,i + 1), 0, 1);
          Vector4f point3(V(0,i + 2), V(1,i + 2), 0, 1);
          Vector4f p1onScreen = model.block(0, model_idx, 4, 4) * point1;
          Vector4f p2onScreen = model.block(0, model_idx, 4, 4) * point2;
          Vector4f p3onScreen = model.block(0, model_idx, 4, 4) * point3;

          animationClicked =  clickOnTriangle(xworld, yworld, p1onScreen[0], p1onScreen[1], p2onScreen[0], p2onScreen[1], p3onScreen[0], p3onScreen[1]);
          if(animationClicked)
          {
            // Set the color of previously selected triangle
            if(vertex_click1 != i){
              if(vertex_click1 > -1){
                C.col(vertex_click1) = originalColor.col(0);
                C.col(vertex_click2) = originalColor.col(1);
                C.col(vertex_click3) = originalColor.col(2);
              }
              // For the selected triangle
              vertex_click1 = i;
              vertex_click2 = i + 1;
              vertex_click3 = i + 2;
              vertexModel = VtoModel(i);

              // initialize animation start vertec offset for V
              Vector4f point1(V(0,i), V(1,i), 0, 1);
              Vector4f p1onScreen = model.block(0, vertexModel, 4, 4) * point1;
              beginX = p1onScreen[0];
              beginY = p1onScreen[1];
              original = translation.block(0, vertexModel, 4,4);

              // Store the color of previously selected triangle
              originalColor.col(0) = C.col(vertex_click1);
              originalColor.col(1) = C.col(vertex_click2);
              originalColor.col(2) = C.col(vertex_click3);

              // Mark selected triangle in green
              C.col(vertex_click1) << 0.0, 1.0, 0.0;
              C.col(vertex_click2) << 0.0, 1.0, 0.0;
              C.col(vertex_click3) << 0.0, 1.0, 0.0;
              VBO_C.update(C);
            }
            break;
          }
        }
      }
      else if (action == GLFW_RELEASE && animationClicked){
        Vector4f point1(V(0,vertex_click1), V(1,vertex_click1), 0, 1);
        Vector4f p1onScreen = model.block(0, vertexModel, 4, 4) * point1;
        endX = p1onScreen[0];
        endY = p1onScreen[1];
        animationClicked = false;
        resetTranslationVariables();
      }
    }
  }

void resetRestBool()
{
  insertionActive = false;
  translationActive = false;
  deletionActive = false;
  colorActive = false;
  animationActive = false;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_RELEASE)
    {
      // Get the size of the window
      int width, height;
      glfwGetWindowSize(window, &width, &height);

      switch (key)
      {
          case GLFW_KEY_I:
              cout << "Insertion mode activated" << endl;
              resetRestBool();
              insertionActive = true;
              resetTranslationVariables();
              break;
          case GLFW_KEY_O:
              cout << "Translation mode activated"  <<endl;
              resetRestBool();
              translationActive = true;
              resetTranslationVariables();
              break;
          case  GLFW_KEY_P:
              cout << "Delete mode activated" << endl;
              resetRestBool();
              deletionActive = true;
              resetTranslationVariables();
              break;
          case  GLFW_KEY_C:
              cout << "Color mode activated" << endl;
              resetRestBool();
              colorActive = true;
              resetTranslationVariables();
              break;
          case  GLFW_KEY_H:
              cout << "Rotate clockwise by 10 degrees"  << endl;
              rotateFun(true);
              break;
          case  GLFW_KEY_J:
              cout << "Rotate counter-clockwise by 10 degrees"  << endl;
              rotateFun(false);
              break;
          case  GLFW_KEY_K:
              cout << "Scale up by 25 percent" << endl;
              scaleFun(true);
              break;
          case  GLFW_KEY_L:
              cout << "Scale down by 25 percent" << endl;
              scaleFun(false);
              break;
          case GLFW_KEY_W:
              cout << "pan up" << endl;
              panFun(1, width,height);
              break;
          case GLFW_KEY_A:
              cout << "pan left" << endl;
              panFun(2, width, height);
              break;
          case GLFW_KEY_S:
              cout << "pan down" << endl;
              panFun(3, width, height);
              break;
          case GLFW_KEY_D:
              cout << "pan right" << endl;
              panFun(4, width, height);
              break;
          case GLFW_KEY_MINUS:
              cout << "Zoom out" << endl;
              zoom(0, width, height);
              break;
          case GLFW_KEY_EQUAL:
              cout << "Zoom in" << endl;
              zoom(1, width, height);
              break;
          case GLFW_KEY_1:
              selectedColor = 1;
              colorVertex();
              break;
          case GLFW_KEY_2:
              selectedColor = 2;
              colorVertex();
              break;
          case GLFW_KEY_3:
              selectedColor = 3;
              colorVertex();
              break;
          case GLFW_KEY_4:
              selectedColor = 4;
              colorVertex();
              break;
          case GLFW_KEY_5:
              selectedColor = 5;
              colorVertex();
              break;
          case GLFW_KEY_6:
              selectedColor = 6;
              colorVertex();
              break;
          case GLFW_KEY_7:
              selectedColor = 7;
              colorVertex();
              break;
          case GLFW_KEY_8:
              selectedColor = 8;
              colorVertex();
              break;
          case GLFW_KEY_9:
              selectedColor = 9;
              colorVertex();
              break;
          case GLFW_KEY_B:
              resetTranslationVariables();
              cout << "Translation in animation mode begins." << endl;
              resetRestBool();
              animationActive = true;
              break;
          case GLFW_KEY_G:
              if(animationActive){
                cout << "Animation activated." << endl;
                translation.block(0, vertexModel, 4,4) = original;
                model.block(0, vertexModel, 4,4) = translation.block(0, vertexModel, 4, 4) * rotation.block(0, vertexModel, 4, 4) * scaling.block(0, vertexModel, 4, 4);
                animate = true;
                t_start = std::chrono::high_resolution_clock::now();
              }
              break;
          default:
              break;
      }
    }
}

int main(void)
{
    GLFWwindow* window;

    // Initialize the library
    if (!glfwInit())
        return -1;

    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Ensure that we get at least a 3.2 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

        // On apple we have to load a core profile with forward compatibility
    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // On apple we have to load a core profile with forward compatibility
    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();


    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    VBO.init();
    V << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    VBO.update(V);

    VBO_C.init();
    C <<
    0, 0, 0,
    0, 0, 0,
    0, 0, 0;
    VBO_C.update(C);

    colors <<
    0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.3, 1.0, 0.0,
    0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.3, 1.0,
    0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 1.0, 0.0, 0.3;

    identityMatrix <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    view <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    model <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    scaling <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    translation <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    rotation <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;


    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar* vertex_shader =
            "#version 150 core\n"
                    "in vec2 position;"
                    "uniform mat4 view;"
                    "uniform mat4 model;"
                    "in vec3 color;"
                    "out vec3 f_color;"
                    "void main()"
                    "{"
                    "    gl_Position = view * model * vec4(position, 0.0, 1.0);"
                    "    f_color = color;"
                    "}";
    const GLchar* fragment_shader =
            "#version 150 core\n"
                    "in vec3 f_color;"
                    "out vec4 outColor;"
                    "uniform vec3 triangleColor;"
                    "void main()"
                    "{"
                    "    outColor = vec4(f_color, 1.0);"
                    "}";

    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    program.init(vertex_shader,fragment_shader,"outColor");
    program.bind();

    // The vertex shader wants the position of the vertices as an input.
    // The following line connects the VBO we defined above with the position "slot"
    // in the vertex shader
    program.bindVertexAttribArray("position",VBO);
    program.bindVertexAttribArray("color",VBO_C);


    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Register the mouse movement callback
    glfwSetCursorPosCallback(window, cursor_pos_callback);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // Bind your VAO (not necessary if you have only one)
        VAO.bind();

        // Bind your program
        program.bind();

        // Example in main_view:
        //Set the uniform value depending on the time difference
        //auto t_now = std::chrono::high_resolution_clock::now();
        //float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
        //glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);

        // Clear the framebuffer
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, view.data());
        // glUniformMatrix4fv Specify the value of a uniform variable 4×4 matrix for the current program object.

        // Insertion display
        if(insertionCount == 1){ // Display 2 lines
          glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, identityMatrix.data());
          glDrawArrays(GL_LINES, (num*3), 2);
        }else if(insertionCount ==  2){ //Display 3 lines
          glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, identityMatrix.data());
          glDrawArrays(GL_LINES, (num*3), 6);
        }else if(insertionCount == 3){
          insertionCount = 0;
        }

        // Draw triangles
        int modelstart = 0;
        for(int tri = 0; tri < (num*3); tri=tri+3){
          if(animate){
            animateTriangle();
          }
          glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, &model(0,modelstart));
          glDrawArrays(GL_TRIANGLES, tri, 3);
          modelstart = modelstart+4;
        }

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();
    VAO.free();
    VBO.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
