/*
Ray Tracer
Julian Biedka
301364984
CMPT 361
*/

#include <GL/glut.h>
#include "primitives.h"

#define RESOLUTION 600		// Screen resolution
#define BOUNCES 16			// Recursion depth
#define PI 3.1415926535

Colour ambient = Colour(0.002, 0.002, 0.005); 	// Background colour
vector<Sphere> Spheres;							// List of spheres
vector<Plane> Planes;							// List of planes
vector<Quadric> Quadrics;						// List of quadrics

double FOV;			// Field of view to be set
double rayMargin;	// Ratio of pixels to world units

// Universal shader function, can be set to anything you want
vector<double> (*ShaderFunc) (vector<double>, Ray, vector<double>, Material, vector<double>, vector<double>);

// Calculate the ratio of pixels to world units
void setFOV(double fov) {
	FOV = fov;
	rayMargin = 2.0 * tan((fov / 2.0) * PI / 180.0) / RESOLUTION;
	cout << "For an fov of " << FOV << " and a resolution of " << RESOLUTION << "\nthere will be " << rayMargin << " units of distance between each ray" << endl;
}

void init() {
	glClearColor(0, 0, 0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, RESOLUTION, 0.0, RESOLUTION);
}

// Prints a vector
void print3D(vector<double> v) {
	cout << v[0] << ", " << v[1] << ", " << v[2];
}

// Returns the smallest positive number
double smallestRoot2D(vector<double> intersections) {
	if (intersections[0] < intersections[1] && intersections[0] >= 0) {
		return intersections[0];
	} else if (intersections[1] < intersections[0] && intersections[1] >= 0) {
		return intersections[1];
	} else {
		return -1;
	}
}

// Calculates a non-negative dot product
double dot3D(vector<double> v1, vector<double> v2) {
	double product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	if (product < 0) {
		return 0;
	} else {
		return product;
	}
}

// Normalizes a vector
vector<double> normalize3D(vector<double> v) {
	double length = sqrt(dot3D(v, v));
	vector<double> norm{ v[0] / length, v[1] / length, v[2] / length };
	return norm;
}

// Returns a 3D vector as the difference of two points
vector<double> fromPoint3D(vector<double> final, vector<double> initial) {
	return {
		final[0] - initial[0],
		final[1] - initial[1],
		final[2] - initial[2]
	};
}

// Calculates a perfectly reflected ray from and incident and normal vector
vector<double> reflection3D(vector<double> incident, vector<double> norm) {
	return normalize3D({
		2 * norm[0] * dot3D(norm, incident) - incident[0],
		2 * norm[1] * dot3D(norm, incident) - incident[1],
		2 * norm[2] * dot3D(norm, incident) - incident[2]
	});	
}

// Calculates a refracted vector based on indicies of refraction
vector<double> refraction3D(Ray ray, vector<double> incident, vector<double> norm, double n1, double n2) {
	double viewDotNorm = dot3D(incident, norm);
	if (n1 > n2 && acos(viewDotNorm) > asin(n2 / n1)) {
		return reflection3D(incident, norm);
	} else {
		double indR = n1 / n2;
		vector<double> incidentRay = normalize3D(ray.GetDirection());

		double cosAngle = viewDotNorm;
		double indexCoefficient = indR * cosAngle - sqrt(1.0 - indR * indR * (1.0 - cosAngle*cosAngle));

		return normalize3D({
			indR * incidentRay[0] + indexCoefficient * norm[0],
			indR * incidentRay[1] + indexCoefficient * norm[1],
			indR * incidentRay[2] + indexCoefficient * norm[2],
		});
	}
}

// Checks to see if an object is in shadow or not, and returns the appropriate shadow value
double castShadow(Ray ray, vector<double> lightPos) {
	double shadowValue = 1;
	
	for (unsigned i = 0; i < Spheres.size(); i++) {
		vector<double> intersections = Spheres[i].Intersect(ray);
		double intersection = smallestRoot2D(intersections);

		// Roundoff correction pretend 0.999 is 1
		if (intersection > 0 && intersection <= 0.999) {
			vector<double> surfacePoint = {
				ray.GetPosition()[0] + intersection * ray.GetDirection()[0],
				ray.GetPosition()[1] + intersection * ray.GetDirection()[1],
				ray.GetPosition()[2] + intersection * ray.GetDirection()[2]
			};

			shadowValue *= Spheres[i].GetMaterial().GetRefr();
		}
	}

	for (unsigned i = 0; i < Quadrics.size(); i++) {
		vector<double> intersections = Quadrics[i].Intersect(ray);
		double intersection = smallestRoot2D(intersections);

		// Roundoff correction pretend 0.999 is 1
		if (intersection > 0 && intersection <= 0.999) {
			vector<double> surfacePoint = {
				ray.GetPosition()[0] + intersection * ray.GetDirection()[0],
				ray.GetPosition()[1] + intersection * ray.GetDirection()[1],
				ray.GetPosition()[2] + intersection * ray.GetDirection()[2]
			};

			shadowValue *= Quadrics[i].GetMaterial().GetRefr();
		}
	}

	for (unsigned i = 0; i < Planes.size(); i++) {
		double intersection = Planes[i].Intersect(ray);

		// Roundoff correction pretend 0.999 is 1
		if (intersection > 0 && intersection <= 0.999) {
			vector<double> surfacePoint = {
				ray.GetPosition()[0] + intersection * ray.GetDirection()[0],
				ray.GetPosition()[1] + intersection * ray.GetDirection()[1],
				ray.GetPosition()[2] + intersection * ray.GetDirection()[2]
			};

			shadowValue *= Planes[i].GetMaterial().GetRefr();
		}
	}

	return shadowValue;
}

// Toon shader variant of the phong lighting equation, see below
vector<double> toonLight(vector<double> colour, Ray ray, vector<double> point, Material surface, vector<double> normal, vector<double> view) {
	Colour dif = surface.GetDif();
	Colour spec = surface.GetSpec();
	
	for (unsigned light = 0; light < Light::GetLights().size(); light++) {
		vector<double> lightPos = Light::GetLights()[light]->GetPosition();

		vector<double> shadowVector = fromPoint3D(point, lightPos);
		Ray shadowRay;
		shadowRay = Ray(lightPos, shadowVector, ray.GetMedium());

		double shadow = castShadow(shadowRay, lightPos);

		Colour lightCol = Light::GetLights()[light]->GetColour();

		vector<double> lightV = normalize3D(fromPoint3D(lightPos, point));

		vector<double> lightRefl = reflection3D(lightV, normal);

		double diffuseDot = dot3D(normal, lightV);
		double specDot = dot3D(lightRefl, view);

		if (diffuseDot < 0.2) {
			diffuseDot = 0.2;
		} else if (diffuseDot < 0.8) {
			diffuseDot = 0.8;
		} else {
			diffuseDot = 1;
		}

		if (specDot < 0.2) {
			specDot = 0;
		} else if (specDot < 0.9) {
			specDot = 0.6;
		} else if (specDot < 0.99) {
			specDot = 0.8;
		} else {
			specDot = 1;
		}

		colour[0] += shadow * lightCol.GetRed() * (
			dif.GetRed() * diffuseDot + spec.GetRed() * pow(specDot, surface.GetExp()));

		colour[1] += shadow * lightCol.GetGreen() * (
			dif.GetGreen() * diffuseDot + spec.GetGreen() * pow(specDot, surface.GetExp()));

		colour[2] += shadow * lightCol.GetBlue() * (
			dif.GetBlue() * diffuseDot + spec.GetBlue() * pow(specDot, surface.GetExp()));
	}

	return colour;
}

// Standard phong lighting equation
vector<double> phongLight(vector<double> colour, Ray ray, vector<double> point, Material surface, vector<double> normal, vector<double> view) {
	Colour dif = surface.GetDif();
	Colour spec = surface.GetSpec();
	
	// Calculate vectors for each light
	for (unsigned light = 0; light < Light::GetLights().size(); light++) {
		vector<double> lightPos = Light::GetLights()[light]->GetPosition();

		// Cast a ray from the light to the object to see if it is in shadow
		vector<double> shadowVector = fromPoint3D(point, lightPos);
		Ray shadowRay;
		shadowRay = Ray(lightPos, shadowVector, ray.GetMedium());

		double shadow = castShadow(shadowRay, lightPos);

		// Get the colour of the light
		Colour lightCol = Light::GetLights()[light]->GetColour();

		// Calculate the light vector
		vector<double> lightV = normalize3D(fromPoint3D(lightPos, point));

		// Calculate the reflected vector
		vector<double> lightRefl = reflection3D(lightV, normal);

		double diffuseDot = dot3D(normal, lightV);
		double specDot = dot3D(lightRefl, view);

		// Add these colours to the colour vector
		colour[0] += shadow * lightCol.GetRed() * (
			dif.GetRed() * diffuseDot + spec.GetRed() * pow(specDot, surface.GetExp()));

		colour[1] += shadow * lightCol.GetGreen() * (
			dif.GetGreen() * diffuseDot + spec.GetGreen() * pow(specDot, surface.GetExp()));

		colour[2] += shadow * lightCol.GetBlue() * (
			dif.GetBlue() * diffuseDot + spec.GetBlue() * pow(specDot, surface.GetExp()));
	}

	return colour;
}

// Recursive raytracing function
Colour trace(Ray ray, unsigned depth) {
	// When we have hit the recursion limit, return the ambient colour
	if (depth == 0)
		return ambient;

	Sphere* s = NULL;
	Plane* p = NULL;
	Quadric* q = NULL;
	double tValue = -1;

	// Check for intersections
	for (unsigned i = 0; i < Spheres.size(); i++) {
		vector<double> intersections = Spheres[i].Intersect(ray);
		double t = smallestRoot2D(intersections);
		
		if ((tValue == -1 && t > 0.0001) || (t > 0.0001 && t < tValue)) {
			tValue = t;
			s = &Spheres[i];
			p = NULL;
			q = NULL;
		}
	}
	
	for (unsigned i = 0; i < Quadrics.size(); i++) {
		vector<double> intersections = Quadrics[i].Intersect(ray);
		double t = smallestRoot2D(intersections);
		
		if ((tValue == -1 && t > 0.0001) || (t > 0.0001 && t < tValue)) {
			tValue = t;
			q = &Quadrics[i];
			p = NULL;
			s = NULL;
		}
	}

	for (unsigned i = 0; i < Planes.size(); i++) {
		double t = Planes[i].Intersect(ray);

		if (t != -1) {
			if ((tValue == -1 && t > 0.0001) || (t > 0.0001 && t < tValue)) {
				tValue = t;
				p = &Planes[i];
				s = NULL;
				q = NULL;
			}
		}
	}

	// If there is a valid intersection, begin the tracing routine
	if (tValue <= 0.0001) {
		return ambient;
	} else if (s != NULL || p != NULL || q != NULL) {
		Material surfMat;
		if (s != NULL) {
			surfMat = s->GetMaterial();
		} else if (p != NULL) {
			surfMat = p->GetMaterial();
		} else {
			surfMat = q->GetMaterial();
		}

		// Get the point on the surface of the intersected object
		vector<double> surfacePoint = {
			ray.GetPosition()[0] + tValue * ray.GetDirection()[0],
			ray.GetPosition()[1] + tValue * ray.GetDirection()[1],
			ray.GetPosition()[2] + tValue * ray.GetDirection()[2]
		};

		// Define the unit vectors
		vector<double> normV, viewV, reflV, refrV;

		// Get the normal vector from the object
		if (s != NULL) {
			normV = normalize3D(s->GetNormal(surfacePoint));
		} else if (p != NULL) {
			normV = normalize3D(p->GetNormal());
		} else {
			normV = normalize3D(q->GetNormal(surfacePoint));
		}

		// Calculate the view vector
		viewV = normalize3D({
			-ray.GetDirection()[0],
			-ray.GetDirection()[1],
			-ray.GetDirection()[2]
			});

		// Calculate reflected vector
		reflV = reflection3D(viewV, normV);

		// Calculate refracted vector
		refrV = refraction3D(ray, viewV, normV, ray.GetMedium(), surfMat.GetIndex());

		// Recurse on the reflected ray
		Ray reflRay = Ray(surfacePoint, reflV, ray.GetMedium());
		Colour reflectedCol = ambient;
		if (surfMat.GetRefl() != 0.0)
			reflectedCol = trace(reflRay, depth - 1);
		
		// Recurse on the refracted ray
		Ray refrRay = Ray(surfacePoint, refrV, surfMat.GetIndex());
		Colour refractedCol = ambient;
		if (surfMat.GetRefr() != 0.0 && dot3D(refrV, refrV) != 0)
			refractedCol = trace(refrRay, depth - 1);

		// Get the surface materials for the ambient calculation
		Colour dif = surfMat.GetDif();
		double reflCoef = surfMat.GetRefl();
		double refrCoef = surfMat.GetRefr();

		vector<double> I{ 0.0, 0.0, 0.0 };

		// Calculate the base ambient and reflected colours
		I = {
			ambient.GetRed() * dif.GetRed() + reflectedCol.GetRed() * reflCoef + refractedCol.GetRed() * refrCoef,
			ambient.GetGreen() * dif.GetGreen() + reflectedCol.GetGreen() * reflCoef + refractedCol.GetGreen() * refrCoef,
			ambient.GetBlue() * dif.GetBlue() + reflectedCol.GetBlue() * reflCoef + refractedCol.GetBlue() * refrCoef
		};

		// Calculate light with the defined shader function
		I = ShaderFunc(I, ray, surfacePoint, surfMat, normV, viewV);

		return Colour(I);
	}
	return ambient;
}

// For every pixel, create and fire a ray
void render() {
	for (int y = 0; y < RESOLUTION; y++) {
		for (int x = 0; x < RESOLUTION; x++) {
			vector<double> position{ 0.0, 0.0, 0.0 };
			vector<double> direction{ rayMargin * ((double)x - ((double)RESOLUTION / 2.0)) , rayMargin * ((double)y - ((double)RESOLUTION / 2.0)), -1 };

			Colour c = trace( Ray(position, direction, 1), BOUNCES);

			glColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
			glBegin(GL_POINTS);

			glVertex2i(x, y);

			glEnd();
		}
	}
	cout << "Done." << endl;
	Light::KillLights(); // Turn out the lights!
}

void sceneOne() {
	cout << "scene 1" << endl;
	setFOV(50);
	ShaderFunc = phongLight;

	new Light(
		Colour(0.9),
		{-2, 3, 0}
	);

	Planes.push_back(Plane(
		0,
		1,
		1,
		15,
		Material(
			Colour(0.2),
			Colour(0.6),
			1,
			0,
			0,
			1
		)
	));

	// Diffuse spheres
	Spheres.push_back(Sphere(
		{-3, 1, -10},
		0.75,
		Material(
			Colour(1, 0, 0),
			Colour(0),
			1,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{-1, 1, -10},
		0.75,
		Material(
			Colour(0, 1, 0),
			Colour(0),
			1,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{1, 1, -10},
		0.75,
		Material(
			Colour(0, 0, 1),
			Colour(0),
			1,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{3, 1, -10},
		0.75,
		Material(
			Colour(1, 1, 1),
			Colour(0),
			1,
			0.0,
			0.0,
			1
		)
	));

	// Specular spheres
	Spheres.push_back(Sphere(
		{-3, -1, -10},
		0.75,
		Material(
			Colour(0.5),
			Colour(0.75),
			1,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{-1, -1, -10},
		0.75,
		Material(
			Colour(0.5),
			Colour(0.75),
			4,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{1, -1, -10},
		0.75,
		Material(
			Colour(0.5),
			Colour(0.75),
			16,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{3, -1, -10},
		0.75,
		Material(
			Colour(0.5),
			Colour(0.75),
			32,
			0.0,
			0.0,
			1
		)
	));
}

void sceneTwo() {
	cout << "scene 2" << endl;
	setFOV(60);
	ShaderFunc = phongLight;

	new Light(
		Colour(0.8, 0.7, 0.5),
		{0, 2, -3}
	);

	new Light(
		Colour(0.9, 0, 0),
		{-5, 5, 0}
	);

	new Light(
		Colour(0, 0, 0.9),
		{5, 5, 0}
	);

	Planes.push_back(Plane(
		0,		// A
		1,		// B
		0,		// C
		0.25,	// D
		Material(
			Colour(0.5),	// Diffuse
			Colour(0.1),	// Specular
			1,				// Specular Exponent
			0.9,			// Reflective
			0,				// Refractive
			1				// Index
		)
	));

	Spheres.push_back(Sphere(
		{-1, 0.5, -5},	// Centre
		0.75,			// Radius
		Material(
			Colour(0, 0, 1),	// Diffuse
			Colour(1),			// Specular
			16,					// Specular exponent
			0.4,				// Reflective
			0.0,				// Refractive
			1					// Index
		)
	));

	Spheres.push_back(Sphere(
		{1, 0.5, -5},
		0.75,
		Material(
			Colour(1, 0, 0),
			Colour(1),
			16,
			0.4,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{ 0, 0.25, -3.5 },
		0.5,
		Material(
			Colour(0, 0, 0),
			Colour(1),
			32,
			1,
			0.0,
			1
		)
	));
}

void sceneThree() {
	cout << "scene 3" << endl;
	setFOV(45);
	ShaderFunc = phongLight;

	new Light(
		Colour(1.5),
		{0, 1, 0}
	);
	
	Planes.push_back(Plane(
		0,		// A
		1,		// B
		0.3,	// C
		4,		// D
		Material(
			Colour(0.65, 0.75, 0.75),	// Diffuse
			Colour(0.1),	// Specular
			1,				// Specular Exponent
			0,				// Reflective
			0,				// Refractive
			1				// Index
		)
	));

	Spheres.push_back(Sphere(
		{-0.75, -0.75, -5},
		0.45,
		Material(
			Colour(0.2, 0, 0),
			Colour(0.9, 0.9, 1),
			64,
			0,
			0.85,
			0.95
		)
	));

	Spheres.push_back(Sphere(
		{0.75, -0.75, -5},
		0.45,
		Material(
			Colour(0, 0.2, 0),
			Colour(0.9, 0.9, 1),
			64,
			0,
			0.85,
			1
		)
	));

	Spheres.push_back(Sphere(
		{-0.75, 0.75, -5},
		0.45,
		Material(
			Colour(0, 0, 0.2),
			Colour(0.9, 0.9, 1),
			64,
			0,
			0.85,
			1.05
		)
	));

	Spheres.push_back(Sphere(
		{0.75, 0.75, -5},
		0.45,
		Material(
			Colour(0.2, 0.2, 0),
			Colour(0.9, 0.9, 1),
			64,
			0,
			0.85,
			1.1
		)
	));

	Spheres.push_back(Sphere(
		{0, 0, -7},
		1.5,
		Material(
			Colour(0, 0.3, 1),
			Colour(0.9, 0.9, 1),
			64,
			0.1,
			0,
			1
		)
	));
}

void sceneFour() {
	cout << "scene 4" << endl;
	setFOV(60);
	ShaderFunc = toonLight;

	new Light(
		Colour(0.8, 1.1, 0.5),
		{5, 6, -2}
	);

	Planes.push_back(Plane(
		0,
		1,
		0,
		1,
		Material(
			Colour(1, 0.5, 0.8),
			Colour(1),
			8,
			0.05,
			0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{-2, 0, -7},			// Centre
		1,						// Radius
		Material(
			Colour(0.3, 0.3, 1),// Diffuse
			Colour(1),			// Specular
			16,					// Specular exponent
			0,				// Reflective
			0,				// Refractive
			1					// Index
		)
	));

	Spheres.push_back(Sphere(
		{-0.5, 1, -11},			// Centre
		2,						// Radius
		Material(
			Colour(1, 0.3, 0.5),// Diffuse
			Colour(1),			// Specular
			4,					// Specular exponent
			0,					// Reflective
			0,					// Refractive
			1					// Index
		)
	));

	Spheres.push_back(Sphere(
		{3, 0.5, -8},			// Centre
		1.5,					// Radius
		Material(
			Colour(0.2),		// Diffuse
			Colour(1),			// Specular
			32,					// Specular exponent
			0.75,				// Reflective
			0,					// Refractive
			1					// Index
		)
	));

	Spheres.push_back(Sphere(
		{0.5, -0.5, -4},			// Centre
		0.5,					// Radius
		Material(
			Colour(0, 0, 0.1),	// Diffuse
			Colour(1),			// Specular
			32,					// Specular exponent
			0,					// Reflective
			0.85,				// Refractive
			0.75					// Index
		)
	));
}

void sceneFive() {
	cout << "scene 5" << endl;
	setFOV(60);
	ShaderFunc = phongLight;

	new Light(
		Colour(0.9, 0.6, 1.1),
		{-10, 2, 2}
	);

	new Light(
		Colour(0.6, 0.9, 1.1),
		{10, 2, 2}
	);
	
	Planes.push_back(Plane(
		0,		// A
		1,		// B
		-0.1,	// C
		5,		// D
		Material(
			Colour(0.5),	// Diffuse
			Colour(0.1),	// Specular
			1,				// Specular Exponent
			0,			// Reflective
			0,				// Refractive
			1				// Index
		)
	));

	Planes.push_back(Plane(
		0,		// A
		-1,		// B
		-0.1,	// C
		5,		// D
		Material(
			Colour(0.5),	// Diffuse
			Colour(0.1),	// Specular
			1,				// Specular Exponent
			0,			// Reflective
			0,				// Refractive
			1				// Index
		)
	));

	Quadrics.push_back(Quadric(
		{ 0, 0, -10 },
		1,
		-1,
		1,
		4,
		Material(
			Colour(0.3),
			Colour(0.5),
			2,
			0.9,
			0.0,
			1
		)
	));

	Quadrics.push_back(Quadric(
		{ -5, 0, -15 },
		1,
		0.25,
		1,
		-1,
		Material(
			Colour(0.8, 0.1, 0.1),
			Colour(0.3),
			2,
			0.0,
			0.0,
			1
		)
	));

	Quadrics.push_back(Quadric(
		{ 0, 0, -15 },
		1,
		0.25,
		1,
		-1,
		Material(
			Colour(0.8, 0.1, 0.1),
			Colour(0.3),
			2,
			0.0,
			0.0,
			1
		)
	));

	Quadrics.push_back(Quadric(
		{ 5, 0, -15 },
		1,
		0.25,
		1,
		-1,
		Material(
			Colour(0.8, 0.1, 0.1),
			Colour(0.3),
			2,
			0.0,
			0.0,
			1
		)
	));

	Spheres.push_back(Sphere(
		{0, 0, -10},
		1.5,
		Material(
			Colour(0, 0, 0),
			Colour(1),
			32,
			0.1,
			1.0,
			0.8
		)
	));

}

void drawRays(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	render();
	glFlush();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(250, 250);
	glutInitWindowSize(RESOLUTION, RESOLUTION);
	glutCreateWindow("Raytracer");

	init();

	if (argc > 1) {
		switch (atoi(argv[1])) {
		case 1: sceneOne(); break;
		case 2: sceneTwo(); break;
		case 3: sceneThree(); break;
		case 4: sceneFour(); break;
		case 5: sceneFive(); break;
		}
	} else {
		cout << "Error: Scene number not specified\n" << "./RayTracer <scene_num>" << endl;
		exit(1);
	}

	drawRays();
	glutMainLoop();
}