#include "light.h"

// Light class that contains a position and colour
// When a new light is created it is automatically added to a static list of lights
vector<Light*> Light::Lights;

Light::Light(Colour col, vector<double> pos) {
	intensity = col;
	position = pos;
	Light::Lights.push_back(this);
}

// Getters
Colour Light::GetColour() const {
	return intensity;
}

vector<double> Light::GetPosition() const {
	return position;
}

vector<Light*> Light::GetLights() {
	return Light::Lights;
}

// Delete all lights in the static list
void Light::KillLights() {
	for (int i = 0; i < Light::Lights.size(); i++) {
		delete Light::Lights[i];
	}
}