#include "ray.h"

// Simple colour class
// Contains rgb values stored as doubles
class Colour {
private:
	double r, g, b;

public:
	Colour();
	Colour(double);
	Colour(double, double, double);
	Colour(vector<double>);

	// Getters and setters
	vector<double> GetColour() const;
	double GetRed() const;
	double GetBlue() const;
	double GetGreen() const;
};