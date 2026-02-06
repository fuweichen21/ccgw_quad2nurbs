#include "CCG_QMS_meshView.h"
#include"../3rdParty/freeglut/include/GL/glut.h"
#include"../3rdParty/MeshLib/core/bmp/RgbImage.h"
#include"../3rdParty/MeshLib/core/viewer/arcball.h"
using namespace CCG_QMSLib;
using namespace MeshLib;

/*------------------ view mesh, texture start-----------------------------*/

/* window width and height */
int win_width, win_height;
int gButton;
int startx, starty;

int shadeFlag = 0;
bool showMesh = true;
bool showUV = false;

/* rotation quaternion and translation vector for the object */
CQrot       ObjRot(0, 0, 1, 0);
CPoint      ObjTrans(0, 0, 0);

/* global mesh */
M mesh;

/* arcball object */
CArcball arcball;

int textureFlag = 2;
/* texture id and image */
GLuint texName;
RgbImage image;
bool hasTexture = false;

//copy frame buffer to an image
/*! save frame buffer to an image "snap_k.bmp"*/
void read_frame_buffer()
{
	static int id = 0;

	GLfloat* buffer = new GLfloat[win_width * win_height * 3];
	assert(buffer);
	glReadBuffer(GL_FRONT_LEFT);
	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_FLOAT, buffer);

	RgbImage image(win_height, win_width);

	for (int i = 0; i < win_height; i++)
		for (int j = 0; j < win_width; j++)
		{
			float r = buffer[(i * win_width + j) * 3 + 0];
			float g = buffer[(i * win_width + j) * 3 + 1];
			float b = buffer[(i * win_width + j) * 3 + 2];

			image.SetRgbPixelf(i, j, r, g, b);
		}
	delete[]buffer;

	char name[256];
	std::ostringstream os(name);
	os << "snape_" << id++ << ".bmp";
	image.WriteBmpFile(os.str().c_str());

}

/*! initialize bitmap image texture */
void initialize_bmp_texture()
{
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &texName);
	glBindTexture(GL_TEXTURE_2D, texName);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,   GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	//	int ImageWidth  = image.GetNumRows();
	//	int ImageHeight = image.GetNumCols();
	int ImageWidth = image.GetNumCols();
	int ImageHeight = image.GetNumRows();
	GLubyte* ptr = (GLubyte*)image.ImageData();

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
		ImageWidth,
		ImageHeight,
		0,
		GL_RGB,
		GL_UNSIGNED_BYTE,
		ptr);

	if (textureFlag == 1)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	else if (textureFlag == 2)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glEnable(GL_TEXTURE_2D);
}

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
	double rot[16];

	glTranslated(ObjTrans[0], ObjTrans[1], ObjTrans[2]);
	ObjRot.convert(rot);
	glMultMatrixd((GLdouble*)rot);
}

/*! the eye is always fixed at world z = +5 */
void setupEye(void) {
	glLoadIdentity();
	gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight()
{
	GLfloat lightOnePosition[4] = { 0, 0, 1, 0 };
	GLfloat lightTwoPosition[4] = { 0, 0, -1, 0 };
	glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
	glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

/*! draw axis */
void draw_axis()
{
	glLineWidth(5.0);
	//x axis
	glColor3f(1.0, 0.0, 0.0);	//red
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(1, 0, 0);
	glEnd();

	//y axis
	glColor3f(0.0, 1.0, 0);		//green
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1, 0);
	glEnd();

	//z axis
	glColor3f(0.0, 0.0, 1.0);	//blue
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 1);
	glEnd();

	glLineWidth(1.0);
}

/*! Draw mesh */
void draw_mesh()
{
	if (hasTexture)
		glBindTexture(GL_TEXTURE_2D, texName);
	glDisable(GL_LIGHTING);
	for (auto pf : It::MFIterator(&mesh))
	{
		glBegin(GL_POLYGON);
		for (auto v : It::FVIterator(&mesh, pf))
		{
			//if (v == NULL)std::cout << "this is wrong!" << std::endl;
			CPoint& pt = v->point();

			/*using halfedges' texture，cfw modify 2024/7/30*/
			H* fvh = mesh.corner(v, pf);
			CPoint2& uv = fvh->uv();

			//CPoint2& uv = v->uv();
			CPoint n;
			switch (shadeFlag)
			{
			case 0:
				n = pf->normal();
				break;
			case 1:
				n = v->normal();
				break;
			}
			glNormal3d(n[0], n[1], n[2]);
			glTexCoord2d(uv[0], uv[1]);

			glColor3f(0.5, 0.5, 0.5);
			glVertex3d(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}
}

/*! draw uv */
void draw_uv()
{
	//glBindTexture(GL_TEXTURE_2D, texName);
	for (auto pf : It::MFIterator(&mesh))
	{
		glBegin(GL_POLYGON);
		for (auto v : It::FVIterator(&mesh, pf))
		{
			CPoint2& uv = v->uv();
			glColor3f(0.0, 0.0, 1.0);
			glVertex3d(uv[0], uv[1], 0.0);
		}
		glEnd();
	}

}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h)
{
	float ar;
	//std::cout << "w:" << w << "\th:" << h << std::endl;
	win_width = w;
	win_height = h;

	ar = (float)(w) / h;
	glViewport(0, 0, w, h);               /* Set Viewport */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// magic imageing commands
	gluPerspective(60.0, /* field of view in degrees */
		ar, /* aspect ratio */
		0.01, /* Z near */
		//-2.5,
		300.0/* Z far */);

	/*if (w <= h)
		glOrtho(-10.0, 10.0, -10.0 * (GLfloat)h / (GLfloat)w, 10.0 * (GLfloat)h /

			(GLfloat)w, -10.0, 10.0);
	else
		glOrtho(-10.0 * (GLfloat)w / (GLfloat)h, 10.0 * (GLfloat)w / (GLfloat)h, -

			10.0, 10.0, -10.0, 10.0);*/


	glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	glutPostRedisplay();

	//win_width = w;
	//win_height = h;
	//glViewport(0, 0, w, h);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//if (w <= h)
	//	glOrtho(-10.0, 10.0, -10.0 * (GLfloat)h / (GLfloat)w, 10.0 * (GLfloat)h /

	//		(GLfloat)w, -10.0, 10.0);
	//else
	//	glOrtho(-10.0 * (GLfloat)w / (GLfloat)h, 10.0 * (GLfloat)w / (GLfloat)h, -

	//		10.0, 10.0, -10.0, 10.0);
	//glMatrixMode(GL_MODELVIEW);
	////glLoadIdentity();
	////glTranslatef(0.0, 0.0, -9.0);


}

/*! display call back function for mesh*/
void display_mesh()
{
	/* clear frame buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//setupLight();
	/* transform from the eye coordinate system to the world system */
	setupEye();
	glPushMatrix();
	/* transform from the world to the ojbect coordinate system */
	setupObject();
	/* draw the mesh */
	if (showMesh)
	{
		draw_mesh();
	}

	if (showUV)
	{
		draw_uv();
	}
	glPopMatrix();
	glutSwapBuffers();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
	printf("w  -  Wireframe Display\n");
	printf("f  -  Flat Shading \n");
	printf("s  -  Smooth Shading\n");
	printf("?  -  Help Information\n");
	printf("esc - quit\n");
}

/*! Keyboard call back function */
void keyBoard_mesh(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '1':
		showMesh = !showMesh;
		break;
	case '2':
		showUV = !showUV;
		break;
	case 'f':
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 0;
		break;
		/*case 'p':
			patchN++;
			break;*/
	case 's':
		//Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 1;
		break;
	case 'w':
		//Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
	case 't':
		textureFlag = (textureFlag + 1) % 3;
		switch (textureFlag)
		{
		case 0:
			glDisable(GL_TEXTURE_2D);
			break;
		case 1:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			break;
		case 2:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			break;
		}
		break;
	case 'o':
		read_frame_buffer();
		break;
	case '?':
		help();
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

/*! setup GL states */
void setupGLstate() {
	GLfloat lightOneColor[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat globalAmb[] = { .1, .1, .1, 1 };
	GLfloat lightOnePosition[] = { .0, 0.0, 1.0, 1.0 };
	GLfloat lightTwoPosition[] = { .0, 0.0, -1.0, 1.0 };

	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
	//glClearColor(0.35, 0.53, 0.70, 0);
	glClearColor(1, 1, 1, 0);
	glShadeModel(GL_SMOOTH);

	//glEnable(GL_LIGHT1);
	//glEnable(GL_LIGHT2);
	//glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);

	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, lightOneColor);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
	glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

/*! mouse click call back function */
void  mouseClick(int button, int state, int x, int y) {
	/* set up an arcball around the Eye's center
	switch y coordinates to right handed system  */

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		gButton = GLUT_LEFT_BUTTON;
		arcball = CArcball(win_width, win_height, x - win_width / 2, win_height - y - win_height / 2);
	}

	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
		startx = x;
		starty = y;
		gButton = GLUT_MIDDLE_BUTTON;
	}

	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
		startx = x;
		starty = y;
		gButton = GLUT_RIGHT_BUTTON;
	}
	return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
	CPoint trans;
	CQrot  rot;

	/* rotation, call arcball */
	if (gButton == GLUT_LEFT_BUTTON)
	{
		rot = arcball.update(x - win_width / 2, win_height - y - win_height / 2);
		ObjRot = rot * ObjRot;
		glutPostRedisplay();
	}

	/*xy translation */
	if (gButton == GLUT_MIDDLE_BUTTON)
	{
		double scale = 10. / win_height;
		trans = CPoint(scale * (x - startx), scale * (starty - y), 0);
		startx = x;
		starty = y;
		ObjTrans = ObjTrans + trans;
		glutPostRedisplay();
	}

	/* zoom in and out */
	if (gButton == GLUT_RIGHT_BUTTON) {
		double scale = 10. / win_height;
		trans = CPoint(0, 0, scale * (starty - y));
		startx = x;
		starty = y;
		ObjTrans = ObjTrans + trans;
		glutPostRedisplay();
	}

}

/*! Normalize mesh
* \param pMesh the input mesh
*/
void normalize_mesh(M* pMesh)
{
	CPoint s(0, 0, 0);
	for (auto v : It::MVIterator(pMesh))
	{
		s = s + v->point();
	}
	s = s / pMesh->numVertices();

	for (auto v : It::MVIterator(pMesh))
	{
		CPoint p = v->point();
		p = p - s;
		v->point() = p;
	}

	double d = 0;
	for (auto v : It::MVIterator(pMesh))
	{
		CPoint p = v->point();
		for (int k = 0; k < 3; k++)
		{
			d = (d > fabs(p[k])) ? d : fabs(p[k]);
		}
	}

	for (auto v : It::MVIterator(pMesh))
	{
		CPoint p = v->point();
		p = p / d;
		v->point() = p;
	}
};

/*! Compute the face normal and vertex normal
* \param pMesh the input mesh
*/
void compute_normal(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		CPoint n(0, 0, 0);
		for (auto pF : It::VCcwFIterator(pMesh, v))
		{
			CPoint p[3];
			CHalfEdge* he = pF->halfedge();
			for (int k = 0; k < 3; k++)
			{
				p[k] = he->target()->point();
				he = he->he_next();
			}

			CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
			pF->normal() = fn / fn.norm();
			n += fn;
		}

		n = n / n.norm();
		v->normal() = n;
	}
};

/*! view mesh (tri or quad)*/
void init_openGL_mesh(int argc, char* argv[])
{
	if (hasTexture)
		image.LoadBmpFile(argv[3]);

	/* glut stuff */
	glutInit(&argc, argv);                /* Initialize GLUT */
	//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Mesh Viewer");	  /* Create window with given title */
	glViewport(0, 0, 800, 600);
	//init();
	glutDisplayFunc(display_mesh);             /* Set-up callback functions */
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyBoard_mesh);
	setupGLstate();

	if (hasTexture)
		initialize_bmp_texture();

	glutMainLoop();                       /* Start GLUT event-processing loop */
}

void CCG_QMSLib::viewMesh(M* pMesh, int argc, char* argv[])
{
	if (argc > 2)
	{
		hasTexture = true;
	}
	mesh = *pMesh;
	normalize_mesh(&mesh);
	compute_normal(&mesh);
	init_openGL_mesh(argc, argv);
}
/*------------------ view mesh, texture end-----------------------------*/

/*------------------ view mesh start-----------------------------*/

M view_triMesh;
M view_quadMesh;
bool showSharpEdges = false;
bool showtripos = false;
bool showtriangle = false;
bool showduu = false;
bool showduv = false;
bool showdvv = false;
bool showTriMesh = true;
bool showQuadMesh = true;

void drawEdges()
{
	glLineWidth(1.0);
	glBegin(GL_LINES);
	for (auto pE : It::MEIterator(&mesh))
	{
		V* p0 = mesh.edgeVertex1(pE);
		V* p1 = mesh.edgeVertex2(pE);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
		glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
	}
	glEnd();
}

void draw_sharp_edges()
{
	glLineWidth(4.0);
	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	for (auto pE : It::MEIterator(&mesh))
	{
		/*if (pE->feature())
		{
			V* p0 = mesh.edgeVertex1(pE);
			V* p1 = mesh.edgeVertex2(pE);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}*/
		if (pE->sharp() || pE->boundary())
		{
			V* p0 = mesh.edgeVertex1(pE);
			V* p1 = mesh.edgeVertex2(pE);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}

	}
	glEnd();
	glLineWidth(1.0);
}

void draw_feature_edges()
{
	glLineWidth(4.0);
	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	for (auto pE : It::MEIterator(&mesh))
	{
		if (pE->feature())
		{
			V* p0 = mesh.edgeVertex1(pE);
			V* p1 = mesh.edgeVertex2(pE);
			glColor3f(0.5f, 0.0f, 1.0f);//purple
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}
		else if (pE->constrinedBoundary())
		{
			V* p0 = mesh.edgeVertex1(pE);
			V* p1 = mesh.edgeVertex2(pE);
			glColor3f(0.0f, 1.0f, 0.0f);//green
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}
	}
	glEnd();
	glLineWidth(1.0);
}

/*对一些特殊点进行可视化*/
void draw_specialVertex()
{
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (auto v : It::MVIterator(&mesh))
	{
		//if (v->id() == 1251 || v->id() == 1252 || v->id() == 1260 || v->id() == 1261 || v->id() == 1262
		//	|| v->id() == 1968 || v->id() == 1969 || v->id() == 1970 || v->id() == 1971 || v->id() == 1976
		//	|| v->id() == 1977 || v->id() == 1978 || v->id() == 1980 || v->id() == 1981 || v->id() == 1982
		//	|| v->id() == 1983 || v->id() == 1986 || v->id() == 1988 || v->id() == 1989)
		//{
		//	if (v->ifSingular())
		//	{
		//		glColor3f(1.0, 0.0, 0.0);//红色
		//	}
		//	else if (v->id() == 1986)
		//	{
		//		glColor3f(1.0, 0.0, 1.0);//品红色（Magenta）
		//	}
		//	else
		//	{
		//		glColor3f(1.0, 1.0, 0.0);//黄色
		//	}
		//	glVertex3f(v->point()[0], v->point()[1], v->point()[2]);
		//}
		if (v->boundary())
		{
			glColor3f(0.5, 0.5, 0.5);//灰色
			glVertex3f(v->point()[0], v->point()[1], v->point()[2]);
		}
		else if (v->ifSharp())
		{
			glColor3f(1.0, 0.0, 1.0);//品红色（Magenta）
			glVertex3f(v->point()[0], v->point()[1], v->point()[2]);
		}
	}
	glEnd();
}

void draw_triMesh()
{
	if (hasTexture)
		glBindTexture(GL_TEXTURE_2D, texName);
	glDisable(GL_LIGHTING);
	for (auto pf : It::MFIterator(&view_triMesh))
	{
		glBegin(GL_POLYGON);
		for (auto v : It::FVIterator(&view_triMesh, pf))
		{
			//if (v == NULL)std::cout << "this is wrong!" << std::endl;
			CPoint& pt = v->point();

			///*使用半边上的参数uv进行贴图，cfw 修改 2024/7/30*/
			//H* fvh = mesh.corner(v, pf);
			//CPoint2& uv = fvh->uv();

			////CPoint2& uv = v->uv();
			//CPoint n;
			//switch (shadeFlag)
			//{
			//case 0:
			//	n = pf->normal();
			//	break;
			//case 1:
			//	n = v->normal();
			//	break;
			//}
			//glNormal3d(n[0], n[1], n[2]);
			//glTexCoord2d(uv[0], uv[1]);

			glColor3f(0.52, 0.52, 0.53);
			glVertex3d(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}

	// Edges
	glLineWidth(1.0);
	glBegin(GL_LINES);
	for (auto pE : It::MEIterator(&view_triMesh))
	{
		V* p0 = view_triMesh.edgeVertex1(pE);
		V* p1 = view_triMesh.edgeVertex2(pE);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
		glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
	}
	glEnd();

	/*画出没有被用来当作采样点的三角网格点*/
	//glPointSize(10.0);
	//glBegin(GL_POINTS);
	//for (auto v : It::MVIterator(&view_triMesh))
	//{
	//	if (!v->ifVisit())
	//	{
	//		glColor3f(1.0, 0.0, 0.0);//红色
	//		glVertex3f(v->point()[0], v->point()[1], v->point()[2]);
	//	}

	//}
	//glEnd();
}
void draw_quadMesh()
{
	if (hasTexture)
		glBindTexture(GL_TEXTURE_2D, texName);
	glDisable(GL_LIGHTING);
	for (auto pf : It::MFIterator(&view_quadMesh))
	{
		if (!pf->tempMark()) continue;
		glBegin(GL_POLYGON);
		for (auto v : It::FVIterator(&view_quadMesh, pf))
		{
			//if (v == NULL)std::cout << "this is wrong!" << std::endl;
			CPoint& pt = v->point();

			///*使用半边上的参数uv进行贴图，cfw 修改 2024/7/30*/
			//H* fvh = mesh.corner(v, pf);
			//CPoint2& uv = fvh->uv();

			////CPoint2& uv = v->uv();
			//CPoint n;
			//switch (shadeFlag)
			//{
			//case 0:
			//	n = pf->normal();
			//	break;
			//case 1:
			//	n = v->normal();
			//	break;
			//}
			//glNormal3d(n[0], n[1], n[2]);
			//glTexCoord2d(uv[0], uv[1]);

			glColor3f(1, 0.0, 0.0);
			glVertex3d(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}
}

void draw_quadMesh_edges()
{
	// Edges
	glLineWidth(5.0);
	glBegin(GL_LINES);
	for (auto pf : It::MFIterator(&view_quadMesh))
	{
		if (!pf->tempMark()) continue;
		for (auto pE : It::FEIterator(&view_quadMesh, pf))
		{
			V* p0 = view_quadMesh.edgeVertex1(pE);
			V* p1 = view_quadMesh.edgeVertex2(pE);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}
	}
	glEnd();
}

/*! display call back function for mesh*/
void display_mesh_triQuad()
{
	/* clear frame buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//setupLight();
	/* transform from the eye coordinate system to the world system */
	setupEye();
	glPushMatrix();
	/* transform from the world to the ojbect coordinate system */
	setupObject();


	/* draw the mesh */
	if (showTriMesh)
	{
		draw_triMesh();
	}
	if (showQuadMesh)
	{
		//draw_quadMesh();
		draw_quadMesh_edges();
	}
	//draw_patch(patchN);
	//draw_patchBoundarys(patchN);

	if (showUV)
	{
		draw_uv();
	}
	glPopMatrix();
	glutSwapBuffers();
}

/*! Keyboard call back function */
void keyBoard_mesh_triQuad(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '1':
		showMesh = !showMesh;
		break;
	case '2':
		showUV = !showUV;
		break;
	case'3':
		showSharpEdges = !showSharpEdges;
		break;
	case'4':
		showTriMesh = !showTriMesh;
		break;
	case'5':
		showQuadMesh = !showQuadMesh;
		break;
	case 'f':
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 0;
		break;
		/*case 'p':
			patchN++;
			break;*/
	case 's':
		//Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 1;
		break;
	case 'w':
		//Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
	case 't':
		textureFlag = (textureFlag + 1) % 3;
		switch (textureFlag)
		{
		case 0:
			glDisable(GL_TEXTURE_2D);
			break;
		case 1:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			break;
		case 2:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			break;
		}
		break;
	case 'o':
		read_frame_buffer();
		break;
	case '?':
		help();
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

/*! view mesh (tri or quad)*/
void init_openGL_mesh_triQuad(int argc, char* argv[])
{
	if (hasTexture)
		image.LoadBmpFile(argv[3]);

	/* glut stuff */
	glutInit(&argc, argv);                /* Initialize GLUT */
	//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Mesh Viewer");	  /* Create window with given title */
	glViewport(0, 0, 800, 600);
	//init();
	glutDisplayFunc(display_mesh_triQuad);             /* Set-up callback functions */
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyBoard_mesh_triQuad);
	setupGLstate();

	if (hasTexture)
		initialize_bmp_texture();

	glutMainLoop();                       /* Start GLUT event-processing loop */
}

void CCG_QMSLib::viewMesh(M* pMesh)
{
	mesh = *pMesh;
	view_quadMesh = *pMesh;
	int argc = 1;
	char** argv = NULL;
	normalize_mesh(&mesh);
	compute_normal(&mesh);
	init_openGL_mesh_triQuad(argc, argv);
}
void CCG_QMSLib::viewMesh(M* triMesh, M* quadMesh, int argc, char* argv[])
{
	if (argc > 3)
	{
		hasTexture = true;
	}
	view_triMesh = *triMesh;
	view_quadMesh = *quadMesh;
	normalize_mesh(&view_triMesh);
	normalize_mesh(&view_quadMesh);
	init_openGL_mesh_triQuad(argc, argv);
}
/*------------------ view mesh end-----------------------------*/

/*------------------ view Scalar On Mesh Face start-----------------------------*/
// scalar field range
static double scalarMin_f, scalarMax_f;

// ---------------------------------------------------------
// auxiliary function: Jet colormap generator
// let t (0.0 ~ 1.0) map to Blue->Cyan->Green->Yellow->Red
// ---------------------------------------------------------
void mapScalarToJetColor(double t, float& r, float& g, float& b) {
	double v = std::max(0.0, std::min(1.0, t));

	// (r,g,b)
	// 0.0 - 0.25: Blue -> Cyan
	// 0.25 - 0.5: Cyan -> Green
	// 0.5 - 0.75: Green -> Yellow
	// 0.75 - 1.0: Yellow -> Red

	if (v < 0.25) {
		r = 0.0;
		g = 4.0 * v;
		b = 1.0;
	}
	else if (v < 0.5) {
		r = 0.0;
		g = 1.0;
		b = 1.0 - 4.0 * (v - 0.25);
	}
	else if (v < 0.75) {
		r = 4.0 * (v - 0.5);
		g = 1.0;
		b = 0.0;
	}
	else {
		r = 1.0;
		g = 1.0 - 4.0 * (v - 0.75);
		b = 0.0;
	}
}

// draw the mesh with face scalar values
void draw_face_scalar()
{
	if (std::abs(scalarMax_f - scalarMin_f) < 1e-12)
		scalarMax_f = scalarMin_f + 1e-12;
	if (hasTexture) glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	for (auto* pf : It::MFIterator(&mesh))
	{
		double e = pf->scalar();
		//normalize to [0,1]
		double t = (e - scalarMin_f) / (scalarMax_f - scalarMin_f);

		// obtain color from jet colormap
		float r, g, b;
		mapScalarToJetColor(t, r, g, b);

		glColor3f(r, g, b);
		glBegin(GL_POLYGON);
		for (auto* v : It::FVIterator(&mesh, pf))
		{
			const CPoint& P = v->point();
			glVertex3d(P[0], P[1], P[2]);
		}
		glEnd();
	}
}

void drawThinEdges()
{
	glLineWidth(1.0f);                    
	glColor3f(0.3f, 0.3f, 0.3f);
	glBegin(GL_LINES);
	for (auto* pE : It::MEIterator(&mesh))
	{
		V* p0 = mesh.edgeVertex1(pE);
		V* p1 = mesh.edgeVertex2(pE);
		glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
		glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
	}
	glEnd();
}

void display_ScalarOnMeshFace()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	setupEye();   // 你的摄像机设置函数
	glPushMatrix();
	setupObject(); // 你的模型变换（旋转/平移）函数

	draw_face_scalar();
	/* draw the mesh */
	if (showMesh)
	{
		drawThinEdges();
	}

	glPopMatrix();
	glutSwapBuffers();
}

/*! Keyboard call back function */
void keyBoard_ScalarOnMeshFace(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '1':
		showMesh = !showMesh;
		break;
	case 'f':
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 0;
		break;
		/*case 'p':
			patchN++;
			break;*/
	case 's':
		//Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 1;
		break;
	case 'w':
		//Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
	case 't':
		textureFlag = (textureFlag + 1) % 3;
		switch (textureFlag)
		{
		case 0:
			glDisable(GL_TEXTURE_2D);
			break;
		case 1:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			break;
		case 2:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			break;
		}
		break;
	case 'o':
		read_frame_buffer();
		break;
	case '?':
		help();
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

void init_openGL_ScalarOnMeshFace(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Vertex Error Visualization");

	glViewport(0, 0, 800, 600);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	// 注册回调
	glutDisplayFunc(display_ScalarOnMeshFace);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyBoard_ScalarOnMeshFace);

	setupGLstate(); // 你的通用 OpenGL 初始化（背景色、灯光等）

	glutMainLoop();
}

void CCG_QMSLib::viewScalarOnMeshFace(M* pMesh, int argc, char* argv[])
{
	mesh = *pMesh;
	// initialize scalar field range
	scalarMin_f = std::numeric_limits<double>::infinity();
	scalarMax_f = -std::numeric_limits<double>::infinity();
	// find the min and max scalar values
	for (auto pf : It::MFIterator(&mesh))
	{
		//pf->scalar() = pf->quality_scaledJocbian();
		pf->scalar() = pf->quality_edgeRatio();
		double val = pf->scalar();
		if (val < scalarMin_f) scalarMin_f = val;
		if (val > scalarMax_f) scalarMax_f = val;
	}
	normalize_mesh(&mesh);
	compute_normal(&mesh);
	init_openGL_ScalarOnMeshFace(argc, argv);
}
/*------------------ view Scalar On Mesh Face end-----------------------------*/

/*------------------ view sampling start -----------------------------*/
std::vector<QB_sampling> total_sts;
int totalPatchNum = 1;
bool showQuadSampling = true;
int fid = 1;//四边形面的id
bool showAllSampling = false;
bool showLocalSampling = false;
bool showBoundarySampling = false;
bool showAllCpts = false;
int showSamplingId = 1;
std::vector<ControlPointData> cpds;
bool showPatchMesh = false;
bool showPatchSampling = false;
int patchId = 1;
int numL = 1;//四边形中采样点的局部id
bool showFeatureEdges = false;

void draw_corners()
{
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (auto v : It::MVIterator(&mesh))
	{
		if (v->ifCorner())
		{
			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(v->point()[0], v->point()[1], v->point()[2]);
		}
	}
	glEnd();
	glPointSize(1.0);
}

/**
 * draw sampling on one quad.
 */
void draw_sampling_one_quad(int fid)
{
	glPointSize(5.0);
	glBegin(GL_POINTS);
	int f_sampling_index = 0;
	for (int index = 0; index < total_sts.size(); index++)
	{
		if (total_sts[index].fId() == fid)
		{
			f_sampling_index = index;
			glColor3f(0.0, 1.0, 1.0);//浅绿色
			glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
		}
	}
	glEnd();
	glPointSize(5.0);
	glLineWidth(5.0);
	F* pf = mesh.idFace(fid);
	glBegin(GL_POLYGON);
	for (auto v : It::FVIterator(&mesh, pf))
	{
		//if (v == NULL)std::cout << "this is wrong!" << std::endl;
		CPoint& pt = v->point();
		CPoint2& uv = v->uv();
		CPoint n;
		switch (shadeFlag)
		{
		case 0:
			n = pf->normal();
			break;
		case 1:
			n = v->normal();
			break;
		}
		//glNormal3d(n[0], n[1], n[2]);
		glNormal3d(1.0, 1.0, 1.0);
		glTexCoord2d(uv[0], uv[1]);

		glColor3f(1.0, 0.5, 0.5);
		glVertex3d(pt[0], pt[1], pt[2]);

	}
	glEnd();
	glPointSize(1.0);
	glLineWidth(1.0);

	///*画当前四边形网格上采样点的参数域*/
	///*画四边形网格面中第一个点*/
	//glPointSize(10.0);
	//glBegin(GL_POINTS);
	//CPoint f_firstvp = mesh.idVertex(total_sts[f_sampling_index].vId())->point();
	//glColor3f(1.0, 0.0, 1.0);//粉色
	//glVertex3d(f_firstvp[0], f_firstvp[1], f_firstvp[2]);
	//glEnd();
	//glPointSize(1.0);


	//F* f = mesh.idFace(fid);
	//H* fh = mesh.faceHalfedge(f);
	//E* fhE = mesh.halfedgeEdge(fh);
	///*指向四边形网格中第一个点的半边应该是fh->nextHalfedge*/
	//H* fhn = mesh.halfedgeNext(fh);
	//E* fhnE = mesh.halfedgeEdge(fhn);
	//CPoint2 cornerA, cornerB, cornerC, cornerD;
	//cornerA[0] = 0.0, cornerA[1] = 0.0;
	//cornerB[0] = fhE->knotInterval(), cornerB[1] = 0.0;
	//cornerC[0] = fhE->knotInterval(), cornerC[1] = fhnE->knotInterval();
	//cornerD[0] = 0.0, cornerD[1] = fhnE->knotInterval();

	///*画当前网格上采样点对应的参数域上的点*/
	//glPointSize(6.0);
	//glBegin(GL_POINTS);
	//for (int index = 0; index < total_sts.size(); index++)
	//{
	//	if (total_sts[index].fId() == fid)
	//	{
	//		QB_sampling sp = total_sts[index];
	//		CPoint2 uv;
	//		uv = cornerA * sp.ws()[0] + cornerB * sp.ws()[1] + cornerC * sp.ws()[2] + cornerD * sp.ws()[3];
	//		glColor3f(1.0, 0.0, 0.0);
	//		glVertex3d(uv[0], uv[1], 0.0);
	//	}
	//}
	//glEnd();
	//glPointSize(1.0);

	//glPointSize(10.0);
	//glBegin(GL_POINTS);
	//glColor3f(1.0, 0.0, 1.0);//粉色
	//glVertex3d(cornerA[0], cornerA[1], 0.0);
	//glEnd();

	//glPointSize(1.0);
	//glColor3f(1.0, 1.0, 0.0);

	//glLineWidth(5.0);
	//glBegin(GL_LINES);
	//{
	//	glColor3f(0.0f, 0.0f, 1.0f);
	//	/*画参数域上的边界线*/
	//	glVertex3f(cornerA[0], cornerA[1], 0.0);
	//	glVertex3f(cornerB[0], cornerB[1], 0.0);
	//	glVertex3f(cornerB[0], cornerB[1], 0.0);
	//	glVertex3f(cornerC[0], cornerC[1], 0.0);
	//	glVertex3f(cornerC[0], cornerC[1], 0.0);
	//	glVertex3f(cornerD[0], cornerD[1], 0.0);
	//	glVertex3f(cornerD[0], cornerD[1], 0.0);
	//	glVertex3f(cornerA[0], cornerA[1], 0.0);
	//	/*画从参数域指向物理域点的连线*/
	//	/*for (int index = 0; index < total_sts.size(); index++)
	//	{
	//		if (total_sts[index].fId() == fid)
	//		{
	//			QB_sampling sp = total_sts[index];
	//			CPoint2 uv;
	//			uv = cornerA * sp.ws()[0] + cornerB * sp.ws()[1] + cornerC * sp.ws()[2] + cornerD * sp.ws()[3];
	//			glColor3f(0.9, 0.9, 0.9);
	//			glVertex3d(uv[0], uv[1], 0.0);
	//			glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
	//		}
	//	}*/

	//}
	//glEnd();
	//glLineWidth(1.0);

}

void draw_sampling(int index)
{
	if (total_sts[index].boundary()) {
		glPointSize(10.0);
		glBegin(GL_POINTS);
		glColor3f(1.0, 0.0, 0.0);//红色
		glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);

		/*---------------------------------------------------------------------------------*/
		/*画四边形网格面中第一个点*/
		CPoint f_firstvp = mesh.idVertex(total_sts[index].vId())->point();
		glColor3f(1.0, 0.0, 1.0);//粉色
		glVertex3d(f_firstvp[0], f_firstvp[1], f_firstvp[2]);
		/*---------------------------------------------------------------------------------*/

		glEnd();
		glPointSize(12.0);
		glBegin(GL_POINTS);
		for (auto id : total_sts[index].index())
		{
			glColor3f(0.0, 1.0, 1.0);//浅绿色
			//glVertex3d(total_cpts[id][0], total_cpts[id][1], total_cpts[id][2]);
			glVertex3f(mesh.idVertex(id)->point()[0], mesh.idVertex(id)->point()[1], mesh.idVertex(id)->point()[2]);
		}
		glEnd();

		glPointSize(1.0);
		glColor3f(1.0, 1.0, 0.0);
		/*glBegin(GL_POLYGON);
		for (auto v : It::FVIterator(&mesh, mesh.idFace(total_sts[index].fId())))
		{
			glVertex3d(v->point()[0], v->point()[1], v->point()[2]);
		}
		glEnd();*/

		/*---------------------------------------------------------------------------------*/
		/*画采样点所在的四边形网格面的边界线*/
		glLineWidth(5.0);
		glBegin(GL_LINES);
		for (auto fh : It::FHEIterator(&mesh, mesh.idFace(total_sts[index].fId())))
		{
			V* p0 = mesh.halfedgeSource(fh);
			V* p1 = mesh.halfedgeTarget(fh);
			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
			glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
		}
		glEnd();
		glLineWidth(1.0);
		//glPointSize(1.0);
		/*---------------------------------------------------------------------------------*/
	}
}

void draw_boundary_sampling()
{
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (int index = 0; index < total_sts.size(); index++)
	{
		if (total_sts[index].boundary())
		{
			//if (!total_sts[index].boundary())continue;
			glColor3f(1.0, 0.0, 0.0);//红色表示待拟合点
			glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
			glColor3f(0.0, 1.0, 0.0);//绿色表示拟合点
			glVertex3d(total_sts[index].pos_approx()[0], total_sts[index].pos_approx()[1], total_sts[index].pos_approx()[2]);
		}
	}
	glEnd();
	glPointSize(1.0);

	glLineWidth(3.0);
	glBegin(GL_LINES);
	{
		for (int index = 0; index < total_sts.size(); index++)
		{
			if (total_sts[index].boundary())
			{
				//if (!total_sts[index].boundary())continue;
				glColor3f(0.0, 0.0, 1.0);//蓝色线连接了红色待拟合点和绿色拟合点
				glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
				glVertex3d(total_sts[index].pos_approx()[0], total_sts[index].pos_approx()[1], total_sts[index].pos_approx()[2]);
			}
		}
	}
	glEnd();
	glLineWidth(1.0);
}

void draw_cpts()
{
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (auto cpt : cpds)
	{
		glColor3f(0.0, 0.0, 1.0); // 蓝色
		glVertex3d(cpt.point1[0], cpt.point1[1], cpt.point1[2]);
		glColor3f(1.0, 0.0, 0.0); // 蓝色
		glVertex3d(cpt.point2[0], cpt.point2[1], cpt.point2[2]);
	}
	glEnd();
	glPointSize(1.0);

	glLineWidth(3.0);
	glBegin(GL_LINES);
	for (auto cpt : cpds)
	{
		glColor3f(0.0, 1.0, 0.0); // 蓝色
		glVertex3d(cpt.point1[0], cpt.point1[1], cpt.point1[2]);
		glVertex3d(cpt.point2[0], cpt.point2[1], cpt.point2[2]);
	}
	glEnd();
}

/**
 * draw patch sampling.
 *
 */
void draw_patchSampling()
{
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (int index = 0; index < total_sts.size(); index++)
	{
		if (mesh.idFace(total_sts[index].fId())->id() == patchId)
		{
			glColor3f(1.0, 0.0, 0.0);
			glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
		}
	}
	glEnd();
	glPointSize(1.0);
}

void draw_sampling()
{
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (int index = 0; index < total_sts.size(); index++)
	{
		//if (!total_sts[index].boundary())continue;
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(total_sts[index].pos()[0], total_sts[index].pos()[1], total_sts[index].pos()[2]);
	}
	glEnd();
	glPointSize(1.0);
}

/**
 * draw patch.
 */
void draw_patchMesh()
{
	if (hasTexture)
		glBindTexture(GL_TEXTURE_2D, texName);

	for (auto pf : It::MFIterator(&mesh))
	{
		if (pf->f_patchId() == patchId)
		{
			glBegin(GL_POLYGON);
			for (auto v : It::FVIterator(&mesh, pf))
			{

				//if (v == NULL)std::cout << "this is wrong!" << std::endl;
				CPoint& pt = v->point();
				CPoint2& uv = v->uv();
				CPoint n;
				switch (shadeFlag)
				{
				case 0:
					n = pf->normal();
					break;
				case 1:
					n = v->normal();
					break;
				}
				//glNormal3d(n[0], n[1], n[2]);
				glNormal3d(1.0, 1.0, 1.0);
				glTexCoord2d(uv[0], uv[1]);

				glColor3f(0.5, 0.5, 0.5);
				glVertex3d(pt[0], pt[1], pt[2]);

			}
			glEnd();
		}
	}
}

void draw_sampling_uv(int index)
{
	QB_sampling sp = total_sts[index];
	F* f = mesh.idFace(sp.fId());
	H* fh = mesh.faceHalfedge(f);
	E* fhE = mesh.halfedgeEdge(fh);
	/*指向四边形网格中第一个点的半边应该是fh->nextHalfedge*/
	H* fhn = mesh.halfedgeNext(fh);
	E* fhnE = mesh.halfedgeEdge(fhn);
	CPoint2 cornerA, cornerB, cornerC, cornerD;
	cornerA[0] = 0.0, cornerA[1] = 0.0;
	cornerB[0] = fhE->knotInterval(), cornerB[1] = 0.0;
	cornerC[0] = fhE->knotInterval(), cornerC[1] = fhnE->knotInterval();
	cornerD[0] = 0.0, cornerD[1] = fhnE->knotInterval();
	CPoint2 uv;
	uv = cornerA * sp.ws()[0] + cornerB * sp.ws()[1] + cornerC * sp.ws()[2] + cornerD * sp.ws()[3];


	glPointSize(10.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(uv[0], uv[1], 0.0);
	glEnd();

	glPointSize(10.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 1.0);
	glVertex3d(cornerA[0], cornerA[1], 0.0);
	glEnd();

	glPointSize(1.0);
	glColor3f(1.0, 1.0, 0.0);

	glLineWidth(5.0);
	glBegin(GL_LINES);
	{
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(cornerA[0], cornerA[1], 0.0);
		glVertex3f(cornerB[0], cornerB[1], 0.0);
		glVertex3f(cornerB[0], cornerB[1], 0.0);
		glVertex3f(cornerC[0], cornerC[1], 0.0);
		glVertex3f(cornerC[0], cornerC[1], 0.0);
		glVertex3f(cornerD[0], cornerD[1], 0.0);
		glVertex3f(cornerD[0], cornerD[1], 0.0);
		glVertex3f(cornerA[0], cornerA[1], 0.0);
	}
	glEnd();
	glLineWidth(1.0);
	//glPointSize(1.0);
}

/*! display call back function for mesh*/
void display_mesh_sampling()
{
	/* clear frame buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//setupLight();
	/* transform from the eye coordinate system to the world system */
	setupEye();
	glPushMatrix();
	/* transform from the world to the ojbect coordinate system */
	setupObject();

	/* draw sharp edges */
	if (showSharpEdges)
	{
		//drawEdges();
		draw_sharp_edges();
		draw_corners();
	}
	/* draw the mesh */
	if (showMesh)
	{
		draw_mesh();
		drawEdges();
		//draw_sampling(showSamplingId);
	}

	if (showQuadSampling)
	{
		//std::cout << "fid" << std::endl;
		draw_sampling_one_quad(fid);
		//draw_sampling_uv_oneQuad(fid);
		/*std::vector<QB_sampling> locl_total_sts;
		for (int index = 0; index < total_sts.size(); index++)
		{
			if (total_sts[index].fId() == fid)
			{
				locl_total_sts.push_back(total_sts[index]);
			}
		}
		if (numL > locl_total_sts.size()) { numL = 0; }
		draw_oneLine_sampling2uv_quad(fid, numL, locl_total_sts);*/
	}

	if (showLocalSampling)
	{
		draw_sampling(showSamplingId);
	}
	if (showBoundarySampling)
	{
		draw_boundary_sampling();
	}

	if (showAllCpts)
	{
		draw_cpts();
	}

	if (showAllSampling)
	{
		draw_sampling();

	}
	//draw_patch(patchN);
	//draw_patchBoundarys(patchN);

	if (showPatchMesh)
	{
		draw_patchMesh();
	}
	if (showPatchSampling)
	{
		draw_patchSampling();
	}

	if (showUV)
	{
		//draw_uv();
		draw_sampling_uv(showSamplingId);
		//draw_sampling_uv();
		//draw_sampling_uv_boundary(showSamplingId);
	}
	glPopMatrix();
	glutSwapBuffers();
}

/*! Keyboard call back function */
void keyBoard_mesh_sampling(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '1':
		showMesh = !showMesh;
		break;
	case '2':
		showUV = !showUV;
		break;
	case'3':
		showSharpEdges = !showSharpEdges;
		break;
	case'4':
		showAllSampling = !showAllSampling;
		break;
	case'5':
		showLocalSampling = !showLocalSampling;
		break;
	case'7':
		showAllCpts = !showAllCpts;
		break;
	case '8':
		showPatchMesh = !showPatchMesh;
		break;
	case '9':
		showPatchSampling = !showPatchSampling;
		break;
	case 'b':
		showBoundarySampling = !showBoundarySampling;
		break;
	case'+':
	{
		fid = (fid + 1) % mesh.numFaces();
		if (fid == 0)
		{
			fid = mesh.numFaces();
		}
		patchId = (patchId + 1) % totalPatchNum;
		if (patchId == 0)
		{
			patchId = totalPatchNum;
		}
		showSamplingId = (showSamplingId + 1) % total_sts.size();
		std::cout << "sampling index: " << showSamplingId << std::endl;
		std::cout << "----uv: " << total_sts[showSamplingId].uv()[0] << " " << total_sts[showSamplingId].uv()[1] << std::endl;
		/*std::cout << "weights: ";
		for (int ii = 0; ii < 4; ii++)
		{
			std::cout << " " << total_sts[showSamplingId].ws()[ii];
		}
		std::cout << std::endl;*/
		break;
	}
	case'-':
	{
		fid = (fid - 1 + mesh.numFaces()) % mesh.numFaces();
		if (fid == 0)
		{
			fid = mesh.numFaces();
		}
		patchId = (patchId - 1 + totalPatchNum) % totalPatchNum;
		if (patchId == 0)
		{
			patchId = totalPatchNum;
		}
		showSamplingId = (showSamplingId - 1 + total_sts.size()) % total_sts.size();
		std::cout << "sampling index: " << showSamplingId << std::endl;
		std::cout << "----uv: " << total_sts[showSamplingId].uv()[0] << " " << total_sts[showSamplingId].uv()[1] << std::endl;
		/*std::cout << "weights: ";
		for (int ii = 0; ii < 4; ii++)
		{
			std::cout << " " << total_sts[showSamplingId].ws()[ii];
		}
		std::cout << std::endl;*/
		break;
	}
	case'A':
	{
		numL = numL - 1;
		break;
	}
	case'a':
	{
		numL = numL + 1;
		break;
	}
	case'6':
		showFeatureEdges = !showFeatureEdges;
		break;
	case'q':
		showQuadSampling = !showQuadSampling;
		break;
	case 'f':
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 0;
		break;
		/*case 'p':
			patchN++;
			break;*/
	case 's':
		//Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		shadeFlag = 1;
		break;
	case 'w':
		//Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
	case 't':
		textureFlag = (textureFlag + 1) % 3;
		switch (textureFlag)
		{
		case 0:
			glDisable(GL_TEXTURE_2D);
			break;
		case 1:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			break;
		case 2:
			glEnable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			break;
		}
		break;
	case 'o':
		read_frame_buffer();
		break;
	case '?':
		help();
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

/*! view mesh (tri or quad)*/
void init_openGL_mesh_sampling(int argc, char* argv[])
{
	if (hasTexture)
		image.LoadBmpFile(argv[2]);

	/* glut stuff */
	glutInit(&argc, argv);                /* Initialize GLUT */
	//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Mesh Viewer");	  /* Create window with given title */
	glViewport(0, 0, 800, 600);
	//init();
	glutDisplayFunc(display_mesh_sampling);             /* Set-up callback functions */
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyBoard_mesh_sampling);
	setupGLstate();

	if (hasTexture)
		initialize_bmp_texture();

	glutMainLoop();                       /* Start GLUT event-processing loop */
}

void CCG_QMSLib::viewSampling(M* pMesh, std::vector<QB_sampling> sts)
{
	mesh = *pMesh;
	total_sts = sts;
	/*std::cout << "sampling index: " << showSamplingId << std::endl;
	std::cout << "weights: ";
	for (int ii = 0; ii < 4; ii++)
	{
		std::cout << " " << total_sts[showSamplingId].ws()[ii] ;
	}
	std::cout<< std::endl;*/

	/*对采样点信息处理一下子*/
	/*for (int i = 0; i < total_sts.size(); i++)
	{
		total_sts[i].fId() = total_sts[i].patchId();
	}*/

	/*获取面片的总个数*/
	for (auto f : It::MFIterator(&mesh))
	{
		if (f->id() > totalPatchNum)
		{
			totalPatchNum = f->id();
		}
	}
	std::cout << " totalPatchNum: " << totalPatchNum << std::endl;
	int argc = 1;
	char** argv = NULL;
	//normalize_mesh_sampling(&mesh);
	//compute_normal(&mesh);
	//normalize_mesh_sampling_uv(&mesh);
	init_openGL_mesh_sampling(argc, argv);
}
/*------------------ view sampling start end-----------------------------*/
