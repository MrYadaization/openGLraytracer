#include "primitives.h"

// Class that contains three type of primitive objects
// Planes, spheres, and generic quadrics

// Basic sphere shape with a position and radius
Sphere::Sphere() {
	position = { 0.0, 0.0, 0.0 };
	radius = 0;
	material = Material();
}

Sphere::Sphere(vector<double> pos, double r, Material mat) {
	position = pos;
	radius = r;
	material = mat;
}

// Check for if a given ray intersects with the object
// Substitute x, y, and z for the rays x(t), y(t), z(t) and solve for t
vector<double> Sphere::Intersect(Ray ray) const {
	double x, y, z, p, q, r, a, b, c, Acoef, Bcoef, Ccoef;
	x = ray.GetPosition()[0];
	y = ray.GetPosition()[1];
	z = ray.GetPosition()[2];

	p = position[0];
	q = position[1];
	r = position[2];

	a = ray.GetDirection()[0];
	b = ray.GetDirection()[1];
	c = ray.GetDirection()[2];

	Acoef = a * a
		+ b * b
		+ c * c;

	Bcoef = 2 * x*a - 2 * p*a
		+ 2 * y*b - 2 * q*b
		+ 2 * z*c - 2 * r*c;

	Ccoef = x * x + p * p - 2 * p*x
		+ y * y + q * q - 2 * q*y
		+ z * z + r * r - 2 * r*z
		- radius * radius;

	double discriminant = Bcoef * Bcoef - 4 * Acoef*Ccoef;

	if (discriminant < 0) {
		return { -1, -1 };
	} else {
		double positiveRoot = (-Bcoef + sqrt(discriminant)) / (2 * Acoef);
		double negativeRoot = (-Bcoef - sqrt(discriminant)) / (2 * Acoef);

		return { positiveRoot, negativeRoot };
	}
}

// Getters
vector<double> Sphere::GetNormal(vector<double> pos) const {
	return {
				pos[0] - position[0],
				pos[1] - position[1],
				pos[2] - position[2]
	};
}

Material Sphere::GetMaterial() const {
	return material;
}

double Sphere::GetRadius() const {
	return radius;
}

vector<double> Sphere::GetPosition() const {
	return position;
}

//// **************** ////

// Quadratic surface of the form A(x - p)^2 + B(y - q)^2 + C(z - r)^2 + D = 0
// Restriced to spheres, ellipsoids, cones, and hyperboloids of one and two sheets
Quadric::Quadric() {
	material = Material();
	A = 0.0;
	B = 0.0;
	C = 0.0;
	D = 0.0;
}

Quadric::Quadric(vector<double> pos, double coefA, double coefB, double coefC, double coefD, Material mat) {
	material = mat;
	p = pos[0];
	q = pos[1];
	r = pos[2];
	
	A = coefA;
	B = coefB;
	C = coefC;
	D = coefD;
}

// Check for if a given ray intersects with the object
// Substitute x, y, and z for the rays x(t), y(t), z(t) and solve for t
vector<double> Quadric::Intersect(Ray ray) const {
	double x, y, z, a, b, c;
	x = ray.GetPosition()[0];
	a = ray.GetDirection()[0];
	
	y = ray.GetPosition()[1];
	b = ray.GetDirection()[1];
	
	z = ray.GetPosition()[2];
	c = ray.GetDirection()[2];

	double alpha, beta, gamma;

	alpha = a*a*A + b*b*B + c*c*C;
	beta = 2*a*A*(x - p) + 2*b*B*(y - q) + 2*c*C*(z - r);
	gamma = A*p*p - 2*A*p*x + A*x*x + 
			B*q*q - 2*B*q*y + B*y*y +
			C*r*r - 2*C*r*z + C*z*z + D;

	double discriminant = beta * beta - 4 * alpha*gamma;

	if (discriminant < 0) {
		return { -1, -1 };
	} else {
		double positiveRoot = (-beta + sqrt(discriminant)) / (2 * alpha);
		double negativeRoot = (-beta - sqrt(discriminant)) / (2 * alpha);

		return { positiveRoot, negativeRoot };
	}
}

// Return the normal vector by calculating the gradient of the shape
// returning the normal to the tangent plane at pos
vector<double> Quadric::GetNormal(vector<double> pos) const {
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
	return { 2*A*(x - p), 2*B*(y - q), 2*C*(z - r) };
}

Material Quadric::GetMaterial() const {
	return material;
}

//// **************** ////

// Basic plane of the form Ax + By + Cz + D = 0
Plane::Plane() {
	material = Material();
	A = 0.0;
	B = 0.0;
	C = 0.0;
	D = 0.0;
}

Plane::Plane(double coefA, double coefB, double coefC, double coefD, Material mat) {
	material = mat;
	A = coefA;
	B = coefB;
	C = coefC;
	D = coefD;
}

// Check for if a given ray intersects with the object
// Substitute x, y, and z for the rays x(t), y(t), z(t) and solve for t
double Plane::Intersect(Ray ray) const {
	double x, y, z, a, b, c;
	x = ray.GetPosition()[0];
	y = ray.GetPosition()[1];
	z = ray.GetPosition()[2];

	a = ray.GetDirection()[0];
	b = ray.GetDirection()[1];
	c = ray.GetDirection()[2];

	double t = -(A*x + B * y + C * z + D) / (A*a + B * b + C * c);

	if (t < 0)
		return -1;
	return t;
}

// Getters
Material Plane::GetMaterial() const {
	return material;
}

vector<double> Plane::GetNormal() const {
	return { A, B, C };
}