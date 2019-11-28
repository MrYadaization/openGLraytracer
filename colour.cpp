#include "colour.h"

// Simple colour class
// Contains rgb values stored as doubles
Colour::Colour() {
	r = 0.0;
	g = 0.0;
	b = 0.0;
}

Colour::Colour(double col) {
	r = col;
	g = col;
	b = col;
}

Colour::Colour(double red, double green, double blue) {
	r = red;
	g = green;
	b = blue;
}

Colour::Colour(vector<double> v) {
	r = v[0];
	g = v[1];
	b = v[2];
}

// Getters and setters
vector<double> Colour::GetColour() const {
	vector<double> c{ r, g, b };
	return c;
}

double Colour::GetRed() const {
	return r;
}

double Colour::GetBlue() const {
	return b;
}

double Colour::GetGreen() const {
	return g;
}