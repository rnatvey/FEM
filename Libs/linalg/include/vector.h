#pragma once
#include <cmath>
class vector2 {

private:

	double X;
	double Y;

	inline void assign(vector2& vec) { this->X = vec.X; this->Y = vec.Y };
	inline void assign(double a) { this->X = a; this->Y = a };
	inline void assign(double x, double y) { this->X = x; this->Y = y };

public:

	//конструкторы
	vector2() : {X = 0.0, Y = 0.0};
	vector2(double a, double b) : {assign(a, b)};
	vector2(double a) : {assign(a)};
	vector2(vector2& vec) : {assign(vec)};

	//функции для работы с векторами

	inline double VxV(const vector2& a, const  vector2& b) { return (a.X * b.X + a.Y * b.Y) } const;
	inline double Norm(const vector2& a) { return std::sqrt(VxV(a,a)) } const;
	inline void UniVector() { this /= Norm(this) };
	inline vector2 UniVector(const vector2& vec) { vector2 res{}; res = vec / Norm(vec); return res } const;

	//**********************************************************************************

	vector2 operator + (const vector2& add) const 
	{
		X + add.X;
		Y + add.Y;
	}

	vector2 operator - (const vector2& sub) const
	{
		X - sub.X;
		Y - sub.Y;
	}

	vector2 operator / (double div) const
	{
		X / div;
		Y / div;
	}

	vector2 operator * (double mul) const
	{
		X * mul;
		Y * mul;
	}

	vector2 operator += (const vector2& add) const
	{
		X += add.X;
		Y += add.Y;
	}

	vector2 operator -= (const vector2& sub) const
	{
		X -= sub.X;
		Y -= sub.Y;
	}

	vector2 operator /= (double div) const
	{
		X /= div;
		Y /= div;
	}

	vector2 operator *= (double mul) const
	{
		X *= mul;
		Y *= mul;
	}

	vector2 operator = (const vector2& vec) const
	{
		X = vec.X;
		Y = vec.Y;
	}

};