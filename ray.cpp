#include "ray.h"

// Simple ray class
// Contains a position, direction, and a refractive index
Ray::Ray() {
	mediumIndex = 1;
    position = {0,0,0};
    direction = {0,0,0};
}

Ray::Ray(vector<double> pos, vector<double> dir, double index) {
	mediumIndex = index;
    position = pos;
    direction = dir;
}

double Ray::GetMedium() const {
	return mediumIndex;
}

// Getters and setters
vector<double> Ray::GetPosition() const { return position; }

vector<double> Ray::GetDirection() const { return direction; }

void Ray::SetPosition(vector<double> p) { position = p; }

void Ray::SetDirection(vector<double> d) { direction = d; }

// Print the values of the ray to the terminal for testing
void Ray::Print() {
	cout << "Pos: (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << ", dir: [" << direction[0] << ", " << direction[1] << ", " << direction[2] << "]" << endl;
}