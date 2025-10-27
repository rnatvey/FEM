#pragma once
#include "BaseElement.h"
// Bate FEM page 339
class PlaneIsoparametric :public BaseElement 
{
private:

	matrix8 stiffness_local;
	matrix8 stiffness_global;
	vector8 displacements_local;
	vector8 displacements_global;

	vector2 point1;
	vector2 point2;
	vector2 point3;
	vector2 point4;

	//matrix3x8 B;
	matrix3 C;
	vector3 deformations;

public:

	virtual auto get_K() const override { return stiffness_global; }
	virtual auto get_f() const override { return dosplacements_global; }

	double coords_interpolated(double s, double r);
	double displacements_interpolated(double s, double r);
	matrix2 Jacobian(double s, double r);
	matrix3x8 B(double s, double r);
	matrix8 F(double s, double r);


};