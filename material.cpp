#include "material.h"

// Material class that contains all relevant variables for the lighting equation
Material::Material() {
	Colour blank = Colour(0.0, 0.0, 0.0);
	diffuse = blank;
	specular = blank;
	reflected = 0.0;
	refracted = 0.0;
	refractionIndex = 1.0;
	specularExp = 1;
}

Material::Material(Colour dif, Colour spec, unsigned exp, double refl, double refr, double refractInd) {
	diffuse = dif;
	specular = spec;
	specularExp = exp;
	reflected = refl;
	refracted = refr;
	refractionIndex = refractInd;
}

// Getters
Colour Material::GetDif() const {
	return diffuse;
}

Colour Material::GetSpec() const {
	return specular;
}

double Material::GetRefl() const {
	return reflected;
}

double Material::GetRefr() const {
	return refracted;
}

unsigned Material::GetExp() const {
	return specularExp;
}

double Material::GetIndex() const {
	return refractionIndex;
}