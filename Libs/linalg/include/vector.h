#pragma once
#include <cmath>
class vector2 {

private:

	double X{};
	double Y{};

	inline void assign(const vector2& vec) { if (this != &vec) { this->X = vec.X; this->Y = vec.Y; } };
	inline void assign(double a) { this->X = a; this->Y = a; };
	inline void assign(double x, double y) { this->X = x; this->Y = y; };

public:

	//������������
	vector2() : X{}, Y{} {}
	vector2(double a, double b) :X{}, Y{} { assign(a, b); }
	vector2(double a) : X{}, Y{} { assign(a); };
	vector2(vector2& vec) : X{}, Y{} { assign(vec); };

	//������� ��� ������ � ���������

	inline double VxV(const vector2& a, const  vector2& b) { return (a.X * b.X + a.Y * b.Y); } const
	inline double Norm(const vector2& a) { return std::sqrt(VxV(a, a)); } const
	inline void UniVector() { assign(*this /= Norm(*this)); };
	inline vector2 UniVector(const vector2& vec) { vector2 res{}; res = vec / Norm(vec); return res; } const

	//**********************************************************************************

	vector2 operator + (const vector2& add) const
	{
		return vector2(this->X + add.X, this->Y + add.Y);
	}

	vector2 operator - (const vector2& sub) const
	{
		return vector2(this->X - sub.X, this->Y - sub.Y);
	}

	vector2 operator / (double div) const
	{
		return vector2(this->X / div, this->Y / div);
	}

	vector2 operator * (double mul) const
	{
		return vector2(this->X * mul, this->Y * mul);
	}

	vector2& operator += (const vector2& add)
	{
		vector2::assign(this->X + add.X, this->Y + add.Y);
		return *this;
	}

	vector2& operator -= (const vector2& sub)
	{
		vector2::assign(this->X - sub.X, this->Y - sub.Y);
		return *this;
	}

	vector2& operator /= (double div)
	{
		vector2::assign(*this / div);
		return *this;
	}

	vector2& operator *= (double mul)
	{
		vector2::assign(*this * mul);
		return *this;
	}

	vector2& operator = (const vector2& vec)
	{
		assign(vec);
		return *this;
	}
};