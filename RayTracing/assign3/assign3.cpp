/*
CSCI 420
Assignment 3 Raytracer

Name: Xiaoting Bi
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>
#include <iostream>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
int col = -1, row = -1;

//the field of view of the camera
#define fov 60.0

#define PI 3.1415926
#define degree2radian(a) (a * PI / 180.0)

#define SAMPLING 1

#define ERR 0.0001

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

typedef struct _Triangle
{
	struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
} Sphere;

typedef struct _Light
{
	double position[3];
	double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

typedef struct _Vector
{
	double x;
	double y;
	double z;
} Vector;

typedef struct _Ray
{
	Vector origin;
	Vector direction;
} Ray;

// image plane
Vector left_down, left_up, right_down, right_up;
// width and height of image plane
double width_image_plane, height_image_plane;

// camera position
Vector camera_pos;

// represents all the pixel on the screen : each pixel contains a color value (r, g, b)
double *** screen;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// function tools of math calculation
// vector normalization	
Vector normalize(Vector v)
{
	double l = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	v.x /= l;
	v.y /= l;
	v.z /= l;

	return v;
}

// dot production of 2 vectors	
double dot_production(Vector a, Vector b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

// cross production of 2 vectors	
Vector cross_production(Vector a, Vector b)
{
	Vector c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;

	return c;
}

// subtraction of 2 vectors: vector a - vector b	
Vector subtraction(Vector a, Vector b)
{
	Vector c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;

	return c;
}

Vector reflection_ray(Vector l, Vector n)
{
	// compute reflection ray: r = 2(l • n) n – l 
	double ln = dot_production(l, n);
	Vector r;
	r.x = 2 * ln * n.x - l.x;
	r.y = 2 * ln * n.y - l.y;
	r.z = 2 * ln * n.z - l.z;

	return r;
}

// Ray position: p(t) = p0 + d t for t > 0 
Vector ray_position_Pt(Ray p0, double t)
{
	Vector p;
	p.x = p0.origin.x + p0.direction.x * t;
	p.y = p0.origin.y + p0.direction.y * t;
	p.z = p0.origin.z + p0.direction.z * t;

	return p;
}

bool ray_sphere_intersection(Ray ray, int num, double &distance, Vector &normal)
{
	double r = spheres[num].radius;
	// Center of sphere = [xc yc zc]T
	// a = 1
	// b = 2 * (xd * (x0-xc) + yd * (y0-yc) + zd * (z0-zc))
	// c = (x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2 - r^2
	double b = 2 * (ray.direction.x * (ray.origin.x - spheres[num].position[0]) + 
		ray.direction.y * (ray.origin.y - spheres[num].position[1]) +
		ray.direction.z * (ray.origin.z - spheres[num].position[2]));
	double c = pow(ray.origin.x-spheres[num].position[0], 2) + 
		pow(ray.origin.y-spheres[num].position[1], 2) +
		pow(ray.origin.z-spheres[num].position[2], 2) - r * r;
	// Solve to obtain t0 and t1: t0,1 = (-b +- sqrt(b^2-4c)) / 2
	double temp = b * b - 4 * c;
	if (temp < 0) // Calculate b2 – 4c, abort if negative 
		return false;
	temp = sqrt(temp);
	double t0 = (-b + temp) / 2;
	double t1 = (-b - temp) / 2;
	if (t0 > 0 && t1 > 0)
		distance = min(t0, t1);
	else if (t0 < 0 && t1 < 0)
		return false;
	else
		distance = max(t0, t1);

	if (distance < ERR) // collision with original sphere; ignore it
		return false;

	// For lighting, calculate unit normal 
	Vector v = ray_position_Pt(ray, distance);
	// 1/r[xi-xc  yi-yc  zi-zc]T
	normal.x = v.x - spheres[num].position[0];
	normal.y = v.y - spheres[num].position[1];
	normal.z = v.z - spheres[num].position[2];
	normal = normalize(normal);

	return true;
}

bool ray_triangle_intersection(Ray ray, int num, double &distance, Vector &normal, double &a, double &b, double &c)
{
	Vector p0, p1, p2; // 3 vertices of a triangle
	p0.x = triangles[num].v[0].position[0];
	p0.y = triangles[num].v[0].position[1];
	p0.z = triangles[num].v[0].position[2];
	p1.x = triangles[num].v[1].position[0];
	p1.y = triangles[num].v[1].position[1];
	p1.z = triangles[num].v[1].position[2];
	p2.x = triangles[num].v[2].position[0];
	p2.y = triangles[num].v[2].position[1];
	p2.z = triangles[num].v[2].position[2];

	// STEP 1: find the intersection of the point with the triangle’s plane
	// normal: n = (p1 − p0) X (p2 − p0)
	normal = cross_production(subtraction(p1, p0), subtraction(p2, p0));
	normal = normalize(normal);
	// if n • d = 0, no intersection, ray parallel to plane
	double nd = dot_production(normal, ray.direction);
	if (nd < 0.0000001 && nd > -0.0000001)
		return false;
	// t = - (o - p) * n / (d n)
	distance = -1 * dot_production(subtraction(ray.origin, p0), normal) / nd;
	// if t <= 0, the intersection is behind ray origin
	if (distance < 0.0000001)
		return false;

	// STEP 2: find out whether that point is inside the triangle or not
	// the signs must follow: 
	// (p1 − p0) X (x − p0) dot n >= 0
	// (p2 − p1) X (x − p1) dot n >= 0
	// (p0 − p2) X (x − p2) dot n >= 0
	// x = o + d t
	Vector x = ray_position_Pt(ray, distance);
	double i = dot_production(cross_production(subtraction(p1, p0), subtraction(x, p0)), normal);
	double j = dot_production(cross_production(subtraction(p2, p1), subtraction(x, p1)), normal);
	double k = dot_production(cross_production(subtraction(p0, p2), subtraction(x, p2)), normal);
	if (i >= 0.0 && j >= 0.0 && k >= 0.0) // inside the triangle
	{
		double area = dot_production(cross_production(subtraction(p1, p0), subtraction(p2, p0)), normal);
		a = i / area;
		b = j / area;
		c = 1.0 - a - b; // k / area;

		return true;
	}

	return false;
}

// compute color w or w/o shadow of the intersection position
void shading(Vector pos, Vector normal, Vector v, Vector &color, Vector kd, Vector ks, double shi)
{
	for (int i = 0; i < num_lights; i ++)
	{
		normal = normalize(normal);
		bool in_shadow = false;
		Vector direction; // shadow ray direction
		direction.x = lights[i].position[0] - pos.x;
		direction.y = lights[i].position[1] - pos.y;
		direction.z = lights[i].position[2] - pos.z;
		direction = normalize(direction);
		Ray shadow_ray;
		shadow_ray.direction = direction;
		shadow_ray.origin = pos;
		double pos2light = sqrt(pow(pos.x-lights[i].position[0], 2) + pow(pos.y-lights[i].position[1], 2)
			+ pow(pos.z-lights[i].position[2], 2));
		// check if there is an intersection with a sphere or triangle
		// check for sphere intersections
		for (int sphere = 0; sphere < num_spheres; sphere ++)
		{
			double dis;
			Vector norm;
			// check if the shadow_ray intersects with any other spheres
			if (ray_sphere_intersection(shadow_ray, sphere, dis, norm))
			{
				Vector pt = ray_position_Pt(shadow_ray, dis);
				double tmp = sqrt(pow(pt.x-pos.x, 2) + pow(pt.y-pos.y, 2) + pow(pt.z-pos.z, 2));
				if (tmp <= pos2light)
					in_shadow = true;
			}
		}

		// check for triangle intersections
		for (int tri = 0; tri < num_triangles; tri ++)
		{
			double dis;
			double a, b, c;
			a = b = c = 0.0;
			Vector norm;
			if (ray_triangle_intersection(shadow_ray, tri, dis, norm, a, b, c))
			{
				Vector pt = ray_position_Pt(shadow_ray, dis);
				double tmp = sqrt(pow(pt.x-pos.x, 2) + pow(pt.y-pos.y, 2) + pow(pt.z-pos.z, 2));
				if (tmp <= pos2light)
					in_shadow = true;
			}

		}

		// no block to the shadow ray, apply a standard Phong model
		if (in_shadow == false) 
		{
			double ln = dot_production(direction, normal);
			if (ln < 0) 
				ln = 0;
			Vector r = normalize(reflection_ray(direction, normal));
			double rv = dot_production(r, v);
			if (rv < 0)
				rv = 0;
			// I = L (kd (l dot n) + ks (r dot v)^a)
			color.x += lights[i].color[0] * (kd.x * ln + ks.x * pow(rv, shi));
			color.y += lights[i].color[1] * (kd.y * ln + ks.y * pow(rv, shi));
			color.z += lights[i].color[2] * (kd.z * ln + ks.z * pow(rv, shi));
		}
	}
}

void ray_tracing(Ray ray, Vector &color)
{
	bool intersect = false;
	double furthest = DBL_MAX;
	ray.direction = normalize(ray.direction);

	// intersection with a sphere
	for (int i = 0; i < num_spheres; i ++)
	{
		double pos2sphere;
		Vector norm;
		if (ray_sphere_intersection(ray, i, pos2sphere, norm))
		{
			if (pos2sphere < furthest)
			{
				norm = normalize(norm);
				intersect = true;
				furthest = pos2sphere;
				Vector p = ray_position_Pt(ray, pos2sphere);
				Vector phone_color;
				phone_color.x = phone_color.y = phone_color.z = 0.0;
				Vector kd, ks;
				kd.x = spheres[i].color_diffuse[0];
				kd.y = spheres[i].color_diffuse[1];
				kd.z = spheres[i].color_diffuse[2];
				ks.x = spheres[i].color_specular[0];
				ks.y = spheres[i].color_specular[1];
				ks.z = spheres[i].color_specular[2];
				double shi = spheres[i].shininess;

				Vector v;
				v.x = -1 * ray.direction.x;
				v.y = -1 * ray.direction.y;
				v.z = -1 * ray.direction.z;
				v = normalize(v);

				shading(p, norm, v, phone_color, kd, ks, shi);
				color.x = phone_color.x * 255.0;
				color.y = phone_color.y * 255.0;
				color.z = phone_color.z * 255.0;
			}			
		}	
	}
	// intersection with a triangle
	for (int i = 0; i < num_triangles; i ++)
	{
		double pos2triangle;
		double a, b, c;
		a = b = c = 0.0;
		Vector norm;
		if (ray_triangle_intersection(ray, i, pos2triangle, norm, a, b, c))
		{
			if (pos2triangle < furthest)
			{
				intersect = true;
				furthest = pos2triangle;
				Vector p = ray_position_Pt(ray, pos2triangle);

				// interpolate normals at position P
				norm.x = triangles[i].v[0].normal[0] * a + triangles[i].v[1].normal[0] * b
					+ triangles[i].v[2].normal[0] * c;
				norm.y = triangles[i].v[0].normal[1] * a + triangles[i].v[1].normal[1] * b
					+ triangles[i].v[2].normal[1] * c;
				norm.z = triangles[i].v[0].normal[2] * a + triangles[i].v[1].normal[2] * b
					+ triangles[i].v[2].normal[2] * c;
				norm = normalize(norm);

				Vector phone_color;
				phone_color.x = phone_color.y = phone_color.z = 0.0;
				Vector kd, ks;
				kd.x = triangles[i].v[0].color_diffuse[0] * a + triangles[i].v[1].color_diffuse[0] * b
					+ triangles[i].v[2].color_diffuse[0] * c;
				kd.y = triangles[i].v[0].color_diffuse[1] * a + triangles[i].v[1].color_diffuse[1] * b
					+ triangles[i].v[2].color_diffuse[1] * c;
				kd.z = triangles[i].v[0].color_diffuse[2] * a + triangles[i].v[1].color_diffuse[2] * b
					+ triangles[i].v[2].color_diffuse[2] * c;
				ks.x = triangles[i].v[0].color_specular[0] * a + triangles[i].v[1].color_specular[0] * b
					+ triangles[i].v[2].color_specular[0] * c;
				ks.y = triangles[i].v[0].color_specular[1] * a + triangles[i].v[1].color_specular[1] * b
					+ triangles[i].v[2].color_specular[1] * c;
				ks.z = triangles[i].v[0].color_specular[2] * a + triangles[i].v[1].color_specular[2] * b
					+ triangles[i].v[2].color_specular[2] * c;
				double shi = triangles[i].v[0].shininess * a + triangles[i].v[1].shininess * b
					+ triangles[i].v[2].shininess * c;

				Vector v;
				v.x = -1 * ray.direction.x;
				v.y = -1 * ray.direction.y;
				v.z = -1 * ray.direction.z;
				v = normalize(v);

				//if (col == 1 && row == 470)
				//{
					//cout << "direction: " << direction.x << "," << direction.y << "," << direction.z << endl;
					//cout << "normal: " << norm.x << "," << norm.y << "," << norm.z << endl;
					//cout << "reflection: " << r.x << "," << r.y << "," << r.z << endl;
				//}

				shading(p, norm, v, phone_color, kd, ks, shi);
				color.x = phone_color.x * 255.0;
				color.y = phone_color.y * 255.0;
				color.z = phone_color.z * 255.0;
			}			
		}	
	}

	if (intersect)
	{
		color.x += ambient_light[0] * 255.0;
		color.y += ambient_light[1] * 255.0;
		color.z += ambient_light[2] * 255.0;
	} else
	{
		color.x = 255.0;
		color.y = 255.0;
		color.z = 255.0;
	}
	color.x = min(max(color.x, 0), 255.0);
	color.y = min(max(color.y, 0), 255.0);
	color.z = min(max(color.z, 0), 255.0);
}

void save_jpg()
{
	Pic *in = NULL;

	in = pic_alloc(640, 480, 3, NULL);
	printf("Saving JPEG file: %s\n", filename);

	memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
	if (jpeg_write(filename, in))
		printf("File saved Successfully\n");
	else
		printf("Error in Saving\n");

	pic_free(in);      

}

//MODIFY THIS FUNCTION
void draw_scene()
{
	unsigned int x, y;
	double yy = left_down.y;
	//simple output
	col = 0;
	for(x = 0; x < HEIGHT * SAMPLING; x ++)
	{
		row = 0;
		double xx = left_down.x;
		for(y = 0; y < WIDTH * SAMPLING; y ++)
		{
			Vector pix;
			pix.x = xx;
			pix.y = yy;
			pix.z = -1;
			Vector direction = subtraction(pix, camera_pos);
			direction = normalize(direction);
			Ray ray2;
			ray2.origin = camera_pos;
			ray2.direction = direction;
			Vector color;
			color.x = color.y = color.z = 0.0;
			ray_tracing(ray2, color);
			screen[x][y][0] = color.x;
			screen[x][y][1] = color.y;
			screen[x][y][2] = color.z;

			xx += width_image_plane / (WIDTH  * SAMPLING);
			row ++;
		}
		yy += height_image_plane / (HEIGHT  * SAMPLING);
		col ++;
	}

	//simple output
	col = 0;
	for(x = 0; x < WIDTH * SAMPLING; x += SAMPLING)
	{
		glPointSize(2.0);  
		glBegin(GL_POINTS);
		row = 0;
		for(y = 0; y < HEIGHT * SAMPLING; y += SAMPLING)
		{
			double r, g, b;
			r = g = b = 0.0;
			// super sampling
			for (int i = 0; i < SAMPLING; i ++)
			{
				for (int j = 0; j < SAMPLING; j ++)
				{
					r += screen[y+j][x+i][0];
					g += screen[y+j][x+i][1];
					b += screen[y+j][x+i][2];
				}
			}
			r = r / pow(SAMPLING, 2.0);
			g = g / pow(SAMPLING, 2.0);
			b = b / pow(SAMPLING, 2.0);

			plot_pixel(col, row, r, g, b);
			row ++;
		}		
		glEnd();
		glFlush();
		col ++;
	}

	if (mode == MODE_JPEG)
	{
		save_jpg();
	}

	printf("Done!\n"); 
	fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
	glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
	glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
	buffer[HEIGHT-y-1][x][0]=r;
	buffer[HEIGHT-y-1][x][1]=g;
	buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
	plot_pixel_display(x,y,r,g,b);
	if(mode == MODE_JPEG)
		plot_pixel_jpeg(x,y,r,g,b);
}

void parse_check(char *expected,char *found)
{
	if(stricmp(expected,found))
	{
		printf("Expected '%s ' found '%s '\n",expected,found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}
}

void parse_doubles(FILE*file, char *check, double p[3])
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
	printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check("rad:",str);
	fscanf(file,"%lf",r);
	printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
	char s[100];
	fscanf(file,"%s",s);
	parse_check("shi:",s);
	fscanf(file,"%lf",shi);
	printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
	FILE *file = fopen(argv,"r");
	int number_of_objects;
	char type[50];
	int i;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file,"%i",&number_of_objects);

	printf("number of objects: %i\n",number_of_objects);

	parse_doubles(file,"amb:",ambient_light);

	for(i=0;i < number_of_objects;i++)
	{
		fscanf(file,"%s\n",type);
		printf("%s\n",type);
		if(stricmp(type,"triangle")==0)
		{

			printf("found triangle\n");
			int j;

			for(j=0;j < 3;j++)
			{
				parse_doubles(file,"pos:",t.v[j].position);
				parse_doubles(file,"nor:",t.v[j].normal);
				parse_doubles(file,"dif:",t.v[j].color_diffuse);
				parse_doubles(file,"spe:",t.v[j].color_specular);
				parse_shi(file,&t.v[j].shininess);
			}

			if(num_triangles == MAX_TRIANGLES)
			{
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if(stricmp(type,"sphere")==0)
		{
			printf("found sphere\n");

			parse_doubles(file,"pos:",s.position);
			parse_rad(file,&s.radius);
			parse_doubles(file,"dif:",s.color_diffuse);
			parse_doubles(file,"spe:",s.color_specular);
			parse_shi(file,&s.shininess);

			if(num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if(stricmp(type,"light")==0)
		{
			printf("found light\n");
			parse_doubles(file,"pos:",l.position);
			parse_doubles(file,"col:",l.color);

			if(num_lights == MAX_LIGHTS)
			{
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else
		{
			printf("unknown type in scene description:\n%s\n",type);
			exit(0);
		}
	}
	return 0;
}

void display()
{

}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0,WIDTH,0,HEIGHT,1,-1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);

	// Camera is at (0,0,0) and points in the negative z-direction 
	camera_pos.x = camera_pos.y = camera_pos.z = 0.0;
	// position of the image plane in frontal view
	double aspect_ratio = (double)WIDTH / HEIGHT;
	double t = tan(degree2radian(fov / 2));
	left_up.x = -aspect_ratio * t;
	left_up.y = t;
	left_down.x = -aspect_ratio * t;
	left_down.y = -t;
	right_up.x = aspect_ratio * t;
	right_up.y = t;
	right_down.x = aspect_ratio * t;
	right_down.y = -t;
	left_up.z = left_down.z = right_up.z = right_down.z = -1;
	// width and height of image plane
	width_image_plane = 2 * aspect_ratio * t;
	height_image_plane = 2 * t;

	// allocate memory for screen pixels
	screen = new double **[HEIGHT * SAMPLING];
	for (int i = 0; i < HEIGHT * SAMPLING; i ++)
	{
		screen[i] = new double *[WIDTH * SAMPLING];
		for (int j = 0; j < WIDTH * SAMPLING; j ++)
			screen[i][j] = new double[3];		
	}
}

void idle()
{
	//hack to make it only draw once
	static int once=0;
	if(!once)
	{
		draw_scene();
		if(mode == MODE_JPEG)
			save_jpg();
	}
	once=1;
}

int main (int argc, char ** argv)
{
	if (argc<2 || argc > 3)
	{  
		printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
		exit(0);
	}
	if(argc == 3)
	{
		mode = MODE_JPEG;
		filename = argv[2];
	}
	else if(argc == 2)
		mode = MODE_DISPLAY;

	glutInit(&argc,argv);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(WIDTH,HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	glutMainLoop();
}
