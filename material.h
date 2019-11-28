#include "light.h"

// Material class that contains all relevant variables for the lighting equation
class Material {
private:
	Colour diffuse, specular;
	double reflected, refracted, refractionIndex;
	unsigned specularExp;

public:
	Material();
	Material(Colour dif, Colour spec, unsigned exp, double refl, double refr, double refractInd);
	
	// Getters
	Colour GetDif() const;
	Colour GetSpec() const;
	double GetRefl() const;
	double GetRefr() const;
	unsigned GetExp() const;
	double GetIndex() const;
};
