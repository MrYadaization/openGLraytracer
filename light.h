#include "colour.h"

// Light class that contains a position and colour
// When a new light is created it is automatically added to a static list of lights
class Light {
private:
	Colour intensity;
	vector<double> position;
	static vector<Light*> Lights;

public:
	Light(Colour col, vector<double> pos);

	// Getters
	Colour GetColour() const;
	vector<double> GetPosition() const;
	static vector<Light*> GetLights();
	
	// Delete all lights in the static list
	static void KillLights();
};