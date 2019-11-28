#pragma once
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Simple ray class
// Contains a position, direction, and a refractive index
class Ray {
private:
	double mediumIndex;
	vector<double> position;
	vector<double> direction;
public:
	Ray();
	Ray(vector<double> pos, vector<double> dir, double index);

	// Getters and setters
	double GetMedium() const;
	vector<double> GetPosition() const;
	vector<double> GetDirection() const;
	void SetPosition(vector<double> p);
	void SetDirection(vector<double> d);

	// Print the values of the ray to the terminal for testing
	void Print();
};