#include "material.h"

// Class that contains three type of primitive objects
// Planes, spheres, and generic quadrics

// Basic sphere shape with a position and radius
class Sphere {
private:
	vector<double> position;
	double radius;
	Material material;

public:
	Sphere();
	Sphere(vector<double>, double, Material);
	
	// Check for if a given ray intersects with the object
	vector<double> Intersect(Ray) const;
	
	// Getters
	vector<double> GetNormal(vector<double>) const;
	Material GetMaterial() const;
	double GetRadius() const;
	vector<double> GetPosition() const;
};

// Quadratic surface of the form A(x - p)^2 + B(y - q)^2 + C(z - r)^2 + D = 0
// Restriced to spheres, ellipsoids, cones, and hyperboloids of one and two sheets
class Quadric {
private:
	double p, q, r, A, B, C, D;
	Material material;

public:
	Quadric();
	Quadric(vector<double>, double, double, double, double, Material);

	// Check for if a given ray intersects with the object
	vector<double> Intersect(Ray) const;

	// Getters
	vector<double> GetNormal(vector<double>) const;
	Material GetMaterial() const;
};

class Plane {
private:
	double A, B, C, D;
	Material material;

// Basic plane of the form Ax + By + Cz + D = 0
public:
	Plane();
	Plane(double, double, double, double, Material);

	// Check for if a given ray intersects with the object
	double Intersect(Ray) const;

	// Getters
	vector<double> GetNormal() const;
	Material GetMaterial() const;
};