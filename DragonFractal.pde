/*
This code demonstrates the basic structure of a recursive method

In this tutorial I extended the basic fractals assigned to include
a much more complicated fractal - the dragon fractal - which implements
more advanced processes and data structures.

The mathematical code which allows for the transformations comes
from a Java linear algebra library I wrote while in highschool. This
library has since been extended and will be available online.

Created by Henry Hammond 2012
 */

int windowWidth = 480;
int windowHeight = 800;


//colors
color screenBackground;

color randomColor() {
  //generate a random color
  return color( (int) random(256), (int) random(256),(int) random(256));
}

void hexagon(float x, float y, float width, float height) {
}

//this is the dragon fractal algorithm 
void dragon(float [][] points, int iteration, int upper) {
  dragon(points,iteration,"R",upper);
}
void dragon(float [][] points, int iteration, String turns,int upperLimit) {

  if(iteration > upperLimit) return;
  if(iteration >= upperLimit) {

    //draw
    stroke(255);
    strokeWeight(1);
    noFill();
    beginShape();
    stroke(255,255,255,255);
    for(int i=0;i<points.length;i++) {
      vertex(points[i][0],points[i][1]);
    }
    endShape(OPEN);

    return;
  }

  int len = points.length-1;
  len = 0;
  float [] top = {

    points[len][0],points[len][1]
  };


  float [][] points2 = new float[points.length*2][2];

  //find center point
  int cx =0;
  int cy =0;
  for(int i=0;i<points.length;i++) {
    cx+= points[i][0];
    cy+= points[i][1];
  }
  cx/=points.length;
  cy/=points.length;
  float [] shift = {
    cx,cy
  };


  //rotate points around center
  points = unshiftPoly(points,shift);
  //points = scalePoly(points,0.5);
  points = rotatePoly(points,-PI/4);
  points = shiftPoly(points,shift);

  float [][] fold;
  float angle = -PI/2;
  if(turns.charAt(turns.length()-1) == 'L') angle*=-1;
  fold = unshiftPoly(points,points[points.length-1]);
  fold = rotatePoly(fold,angle);
  fold = shiftPoly(fold,points[points.length-1]);

  for(int i=0,upper=points.length;i<points.length;i++) {
    points2[i] = points[i];
    points2[i+upper]=fold[upper-i-1];
  }


  float scaler = 1/1.4142356;

  points2 = unshiftPoly(points2,shift);
  points2 = scalePoly(points2,scaler);
  points2 = shiftPoly(points2,shift);

  //find center of fractal
  cx =0;
  cy =0;
  for(int i=0;i<points2.length;i++) {
    cx+= points2[i][0];
    cy+= points2[i][1];
  }
  cx/=points2.length;
  cy/=points2.length;


  float [] d = new float[2];
  d[0]=0;
  d[1]=cy-height/2;
  if(d[1] > 0) {
    points2 = unshiftPoly(points2,d);
  }
  else {
    points2 = unshiftPoly(points2,d);
  }
  d[0]=cx-width/2;
  d[1]=0;
  
  
  if(d[0] > 0) {
    points2 = unshiftPoly(points2,d);
  }
  else {
    points2 = unshiftPoly(points2,d);
  }

  String s = "";
  if(turns.charAt(turns.length()-1)=='R') {
    s = turns+"R";
  }
  else {
    s = turns+"L";
  }


  dragon(points2,iteration+1,s,upperLimit);
  //dragon(points);
}




float[] shiftPoint(float[] originalPoint, float[] amount) {
  float[] n = { 
    originalPoint[0] + amount[0], originalPoint[1] + amount[1]
  };

  return n;
}

float[] unshiftPoint(float[] originalPoint, float[] amount) {
  float[] n = { 
    originalPoint[0] - amount[0], originalPoint[1] - amount[1]
  };

  return n;
}

float[] rotatePoint(float[] point, float angle) {
  float     a = angle(point);
  float[][] r = rotation(angle);

  return multiply(point, r);
}

float[] scalePoint(float[] point, float magnitude) {
  return multiply(point, ScaleMatrix(magnitude));
}

float[] stretchPoint(float[] point, float x, float y) {
  return multiply(point, stretch(x, y));
}

float[][] rotatePoly(float[][] poly, float angle) {
  float[][] r = rotation(angle);
  float[][] n = new float[poly.length][poly[0].length];

  for (int i = 0; i < n.length; i++) {
    n[i] = multiply(poly[i], r);
  }

  return n;
}

float[][] scalePoly(float[][] poly, float amount) {
  return stretchPoly(poly, amount, amount);
}

float[][] stretchPoly(float[][] poly, float x, float y) {
  float[][] n = new float[poly.length][poly[0].length];

  for (int i = 0; i < n.length; i++) {
    n[i] = stretchPoint(poly[i], x, y);
  }

  return n;
}

float[][] shiftPoly(float[][] poly, float[] amount) {
  float[][] n = new float[poly.length][poly[0].length];

  for (int i = 0; i < n.length; i++) {
    n[i] = shiftPoint(poly[i], amount);
  }

  return n;
}

float[][] unshiftPoly(float[][] poly, float[] amount) {
  float[][] n = new float[poly.length][poly[0].length];

  for (int i = 0; i < n.length; i++) {
    n[i] = unshiftPoint(poly[i], amount);
  }

  return n;
}

float angleFromP2P(float[] p1, float[] p2) {
  float[] n = subtract(p1, p2);

  return angle(n);
}

float [] moveByAngle(float angle, float distance) {
  float [] n = new float [2];
  n[0] = distance*cos(angle);
  n[1] = distance*sin(angle);

  return n;
}

float projectOnX(float[] vector) {
  return vector[0];
}
float projectOnX(float magnitude,float angle) {
  return magnitude*cos(angle);
}
float projectOnY(float magnitude,float angle) {
  return magnitude*sin(angle);
}

/**
 *
 * @param col column vector to be converted
 * @return converted column vector in form of {a...n}
 */
float[] toRowVector(float[][] col) {
  assert col[0].length == 1;

  float[] n = new float[col.length];

  for (int i = 0; i < col.length; i++) {
    n[i] = col[i][0];
  }

  return n;
}

/**
 *
 * Convert row vector to a column vector
 *
 * @param row vector to be converted
 * @return column vector in from of {a},...,{n}
 */
float[][] toColVector(float[] row) {
  float[][] n = new float[row.length][1];

  for (int i = 0; i < row.length; i++) {
    n[i][0] = row[i];
  }

  return n;
}

/**
 * Convert a row vector into a standardized form to allow easier calculations
 * @param row vecttor to be converted
 * @return formatted row vector
 */
float[][] stdRowVector(float[] row) {
  return transpose(toColVector(row));
}

/**
 * Calculate dot product of two vectors
 * @param vectorA input vector 1
 * @param vectorB input vector 2
 * @return numerical calculation of dot product
 */
float dotProduct(float[] vectorA, float[] vectorB) {
  return vectorA[0] * vectorB[0] + vectorA[1] * vectorB[1];
}

/**
 * Calculate angle between two vectors.
 * <br /><b>Note: Will only return values between 0 and PI</b>
 * @param vectorA
 * @param vectorB
 * @return angle between 0 and PI
 */
float angleBetween(float[] vectorA, float[] vectorB) {
  float dot = dotProduct(vectorA, vectorB);

  if (dot == 0) {    // tan pi/2 = undefined...
    return PI / 2;
  }

  return acos(dot / (magnitude(vectorA) * magnitude(vectorB)));
}

/**
 * Calculate angle between x axis and a vector or point
 * @param vector
 * @return angle between 0 and 2PI
 */
float angle(float[] vector) {
  float[] baseline = {
    1, 0
  };
  float a = angleBetween(vector, baseline);

  if (vector[1] <= baseline[1]) {
    a *= -1;
  }

  return a;
}

/**
 * Length of a vector
 * @param vector
 * @return numeric value of a vector's length
 */
float magnitude(float[] vector) {
  float v = 0;

  for (int i = 0; i < vector.length; i++) {
    v += vector[i] * vector[i];
  }

  return sqrt(v);
}

/**
 * Multiply two row vectors in non formatted form.
 * <br />This is just a bit faster than using formatted form
 * @param v1
 * @param v2
 * @return
 */
float[] multiply(float[] v1, float v2[]) {
  assert v1.length == v2.length;

  float[] n = new float[1];

  for (int i = 0; i < n.length; i++) {
    n[i] = 0;

    for (int j = 0; j < v1.length; j++) {
      n[i] += v1[j] * v2[j];
    }
  }

  return n;
}

/**
 * Multiply any two matrices.
 * @param m1
 * @param m2
 * @return
 */
float[][] multiply(float[][] m1, float[][] m2) {
  if ((m1[0].length == m2.length) != true) {
    return null;
  }

  float[][] n = new float[m1.length][m2[0].length];

  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      n[r][c] = 0;
      for (int i = 0; i < m1[0].length; i++) {
        n[r][c] += m1[r][i] * m2[i][c];
      }
    }
  }

  return n;
}

/**
 * Multiply vector by a matrix
 * @param v
 * @param m
 * @return
 */
float[] multiply(float[] v, float[][] m) {
  assert v.length == m[0].length;

  float[] n = new float[v.length];

  for (int i = 0; i < n.length; i++) {
    n[i] = 0;

    for (int j = 0; j < m.length; j++) {
      n[i] += v[j] * m[i][j];
    }
  }

  return n;
}

/**
 * Multiply a matrix by a scalar
 * @param m any matrix
 * @param x any scalar
 * @return
 */
float[][] multiply(float[][] m, float x) {
  float[][] n = m;

  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      n[r][c] *= x;
    }
  }

  return m;
}

float[] multiply(float[] m, float x) {
  float[] n = new float[m.length];
  for (int i = 0; i < m.length; i++) {
    n[i] = m[i] * x;
  }

  return n;
}

/**
 * Generate a rotation matrix
 * @param angle of rotation
 * @return a new rotation matrix
 */
float[][] rotation(float angle) {
  float cos = cos(angle);
  float sin = sin(angle);
  float[][] n = {
    {
      cos, -sin
    }
    , {
      sin, cos
    }
    ,
  };

  return n;
}

/**
 * ScaleMatrix a matrix to a larger size
 * @param magnitude
 * @return a new ScaleMatrix matrix
 */
float[][] ScaleMatrix(float magnitude) {
  return stretch(magnitude, magnitude);
}

float[][] stretch(float x, float y) {
  float[][] n = {
    {
      x, 0
    }
    ,
    {
      0, y
    }
  };
  return n;
}

/**
 * Shear matrix on x axis
 * @param shear
 * @return a new shear matrix
 */
float[][] ShearYMatrix(float shear) {
  float[][] n = {
    {
      1, 0
    }
    , {
      shear, 1
    }
    ,
  };

  return n;
}

/**
 * Shear matrix on y axis
 * @param shear
 * @return a new shear matrix
 */
float[][] ShearXMatrix(float shear) {
  float[][] n = {
    {
      1, shear
    }
    , {
      0, 1
    }
    ,
  };

  return n;
}

/**
 * flip a matrix in the x direction
 * @return new transform matrix
 */
float[][] flipX() {
  float[][] n = {
    {
      -1, 0
    }
    , {
      0, 1
    }
    ,
  };

  return n;
}

/**
 * flip a matrix in the y direction
 * @return new transform matrix
 */
float[][] flipY() {
  float[][] n = {
    {
      1, 0
    }
    , {
      0, -1
    }
    ,
  };

  return n;
}

/**
 * Find and return the determinant of a matrix
 * @param m matrix
 * @return determinant
 */
float det(float[][] m) {
  //trace(m);
  assert m.length == m[0].length;

  if (m.length == 2) {
    float n = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    return n;
  } 
  else {
    float n = 0;
    float tmp = 0;

    for (int col = 0; col < m.length; col++) {
      tmp = 1;

      for (int h = 0; h < m.length; h++) {
        tmp *= m[h][(col + h) % m.length];
      }

      n += tmp;
    }

    for (int col = 0; col < m.length; col++) {
      tmp = 1;

      for (int h = 0; h < m.length; h++) {
        tmp *= m[h][(m.length + (col - h)) % m.length];
      }

      n -= tmp;
    }

    return n;
  }
}

float[][] invert(float[][] m) {
  assert isSquare(m);

  if (det(m) != 0) {

    float[][] n = multiply(adj(m), 1 / det(m));
    return n;
  }
  return null;
  //return n;
}

float[][] adj(float[][] m) {

  float[][] n = transpose(m);
  int[] coords = new int[2];
  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n.length; c++) {
      coords[0] = c;
      coords[1] = r;
      n[r][c] = det(minor(coords, m));
    }
  }

  int i = -1;
  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      i *= -1;
      n[r][c] *= i;
    }
  }
  return n;
}

float[][] minor(int[] coords, float[][] m) {
  float[][] minor = new float[2][2];
  int rows = 0;
  int cols = 0;

  for (int r = 0; r < m.length; r++) {
    if (r != coords[0]) {
      for (int c = 0; c < m[0].length; c++) {
        if (c != coords[1]) {
          minor[rows][cols] = m[r][c];
          cols++;
          if (cols >= 2) {
            cols = 0;
            rows++;
          }
        }
      }
    }
  }

  return minor;
}

/**
 * Reflect data in an axis
 * @param angle angle of line of reflection
 * @return
 */
float[][] reflection(float angle) {

  float sin = sin(2 * angle);
  float cos = cos(2 * angle);

  float[][] n = {
    {
      cos, sin
    }
    , {
      sin, -cos
    }
  };

  return n;
}

/**
 * addition two matrices
 * @param a
 * @param b
 * @return Summative matrix of A and B
 */
float[][] addition(float[][] a, float[][] b) {
  assert (a.length == b.length) && (a[0].length == b[0].length);

  float[][] n = a;

  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      n[r][c] += b[r][c];
    }
  }

  return n;
}

/**
 * addition any two vectors or points
 * @param a
 * @param b
 * @return The new point
 */
float[] addition(float[] a, float[] b) {
  assert a.length == b.length;

  float[] n = a.clone();

  for (int i = 0; i < n.length; i++) {
    n[i] += b[i];
  }

  return n;
}

/**
 * Subtract vector a from vector b
 * @param a
 * @param b
 * @return new vector (a-b)
 */
float[] subtract(float[] a, float[] b) {
  assert a.length == b.length;

  float[] n = a.clone();

  for (int i = 0; i < n.length; i++) {
    n[i] -= b[i];
  }

  return n;
}

/**
 * Subtract matrix a from matrix b
 * @param a
 * @param b
 * @return new matrix
 */
float[][] subtract(float[][] a, float[][] b) {
  assert (a.length == b.length) && (a[0].length == b[0].length);

  float[][] n = a;

  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      n[r][c] -= b[r][c];
    }
  }

  return n;
}

/**
 * get width of a matrix
 * @param m matrix
 * @return
 */
int width(float[][] m) {
  return m[0].length;
}

/**
 * get height of a matrix
 * @param m matrix
 * @return
 */
int height(float[][] m) {
  return m.length;
}

/**
 * Check if a matrix is square or not.
 * @param m
 * @return
 */
boolean isSquare(float[][] m) {
  return (m.length == m[0].length);
}

/**
 * Generate transpose of any row vector m
 * @param m
 * @return
 */
float[][] transpose(float[] m) {
  return toColVector(m);
}

/**
 * Generate transpose matrix of any input matrix
 * @param m
 * @return
 */
float[][] transpose(float[][] m) {
  float[][] n = new float[m[0].length][m.length];

  for (int r = 0; r < n.length; r++) {
    for (int c = 0; c < n[0].length; c++) {
      n[r][c] = m[c][r];
    }
  }

  return n;
}

/**
 * Check if two matrices are equal
 * @param m1
 * @param m2
 * @return
 */
boolean isEqual(float[][] m1, float[][] m2) {
  if ((m1.length != m2.length) || (m1[0].length != m2[0].length)) {
    return false;
  }

  for (int r = 0; r < m1.length; r++) {
    for (int c = 0; c < m1[0].length; c++) {
      if (m1[r][c] != m2[r][c]) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Check if any matrix is symmetrical
 * @param m
 * @return
 */
boolean symmetric(float[][] m) {
  return isEqual(m, transpose(m));
}

static float distance(float[] a, float[] b) {

  float d = sqrt(pow(b[0] - a[0], 2) + pow(b[1] - a[1], 2));

  return d;
}

void sierpinskiTriangle(float x, float y, float width, float height) {

  if(width < 5) return;

  color fillColor = randomColor();
  fill(fillColor);
  //noStroke();
  triangle(x,y+height,x+width,y+height,x+width/2,y);

  sierpinskiTriangle(x,y+height/2, width/2, height/2);
  sierpinskiTriangle(x+width/2,y+height/2, width/2, height/2);
  sierpinskiTriangle(x+width/4,y, width/2, height/2);
}

void sierpinskiCarpet(int x, int y, int width, int height) {

  if(width < 5) return;

  color fillColor = randomColor();
  stroke(#000000);
  strokeWeight(1);
  fill(fillColor);

  fillColor = randomColor();

  rect(x+width/3,y+height/3,width/3,height/3);

  int r = 0;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
}
void sierpinskiCarpet(int x, int y, int width, int height, color c) {

  if(width < 5) return;

  color fillColor = randomColor();
  stroke(#000000);
  strokeWeight(1);
  fill(c);
  rect(x+width/3,y+height/3,width/3,height/3);

  int r = 0;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
  sierpinskiCarpet(x + 0*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 1*width/3, y + r*height/3, width/3, height/3,fillColor);
  sierpinskiCarpet(x + 2*width/3, y + r*height/3, width/3, height/3,fillColor);
  r++;
}


void hawaiianEarings(int x, int y, int width, int height) {

  //basis case
  if(width < 5) return;
  color fillColor = randomColor(); 
  ellipseMode(RADIUS); //use circle center and radius
  //This ellseMode means in ellipse(x,y,w,h) x and y refer to center of the 
  //shape, not the upper left corner
  //and w and h refer to radius not diameter

  stroke(#000000);
  strokeWeight(1);
  fill(fillColor);
  //draw circular grid lines
  ellipse(x, y, width/2, height/2); 

  //recursive cases   
  hawaiianEarings(x-width/4,y,width/2,height/2);
  hawaiianEarings(x+width/4,y,width/2,height/2);
  hawaiianEarings(x,y-width/4,width/2,height/2);
  hawaiianEarings(x,y+width/4,width/2,height/2);
}

void setup() {

//  /orientation(PORTRAIT);

  size(displayWidth, displayHeight, P3D);
  windowWidth = displayWidth;
  windowHeight = displayHeight;
  //size(windowWidth,windowHeight);
  //smooth();
  //initialize colors
   
  frameRate(1); //frames per second
  screenBackground = color(0,0,0);
  background(screenBackground);
  hexagon(0,0,windowWidth,windowHeight);
  //sierpinskiTriangle(0,0,windowWidth,windowHeight);
  //sierpinskiCarpet(0,0,windowWidth,windowHeight);
  //hawaiianEarings(windowWidth/2, windowHeight/2, windowWidth, windowHeight);
  
  float [][] points = {
    {width/4,height/2},
    //{width/2+500,height/2},
    //{width/2+180,height/2+100},
    //{width/2,height/2+900},
    //{width/2,height/2},
    {width*3/4,height/2},
    
  };
  
  //dragon(points,0,"R",12);
  
}

int current = 0;

void draw() {
  
  
  
  noStroke();
  fill(0,0,0,40);
  //rect(0,0,width,height);
  float [][] points = {
    {width/4,height/2},
    //{width/2+500,height/2},
    //{width/2+180,height/2+100},
    //{width/2,height/2+900},
    //{width/2,height/2},
    {width*3/4,height/2},
    
  };
  /*
  if(current <= 14){
    background(screenBackground);
    dragon(points,0,"R",current+1);
    current++;
  }
  */
  //intentionally empty
}

