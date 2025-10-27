#pragma once
#include <cmath>
class vector2 {

private:

	inline void assign(const vector2& vec) { if (this != &vec) { this->X = vec.X; this->Y = vec.Y; } };
	inline void assign(double a) { this->X = a; this->Y = a; };
	inline void assign(double x, double y) { this->X = x; this->Y = y; };

public:

	double X{ 0.0 };
	double Y{ 0.0 };

	//**********************************************************************************
	vector2() : X{}, Y{} {}
	vector2(double a, double b) :X{}, Y{} { assign(a, b); }
	vector2(double a) : X{}, Y{} { assign(a); };
	vector2(vector2& vec) : X{}, Y{} { assign(vec); };

	//**********************************************************************************

	inline double VxV(const vector2& a, const  vector2& b) { return (a.X * b.X + a.Y * b.Y); } const
	inline double Norm(const vector2& a) { return std::sqrt(VxV(a, a)); } const
	inline void UniVector() { assign(*this /= Norm(*this)); };
	inline const vector2& UniVector(const vector2& vec) { vector2 res{}; res = vec / Norm(vec); return res; } const
	inline double get_x() const { return this->X; } 
	inline double get_y() const { return this->Y; } 
	inline vector2 get_vector() {};
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

//**********************************************************************************
//**********************************************************************************

class vector8 {

public:

	double X_i{ 0.0 };
	double Y_i{ 0.0 };
	double X_j{ 0.0 };
	double Y_j{ 0.0 };
	double X_k{ 0.0 };
	double Y_k{ 0.0 };
	double X_l{ 0.0 };
	double Y_l{ 0.0 };

private:

	inline void assign(const vector8& vec) { if (this != &vec) {
		this->X_i = vec.X_i; this->Y_i = vec.Y_i;
		this->X_j = vec.X_j; this->Y_j = vec.Y_j;
		this->X_k = vec.X_k; this->Y_k = vec.Y_k;
		this->X_l = vec.X_l; this->Y_l = vec.Y_l;
	} };
	inline void assign(double a) { 
		this->X_i = a; this->Y_i = a;
		this->X_j = a; this->Y_j = a;
		this->X_k = a; this->Y_k = a;
		this->X_l = a; this->Y_l = a;
	};
	inline void assign(double xi, double yi, double xj, double yj, double xk, double yk, double xl, double yl) {
		this->X_i = xi; this->Y_i = yi;
		this->X_j = xj; this->Y_j = yj;
		this->X_k = xk; this->Y_k = xk;
		this->X_l = xl; this->Y_l = xl;
	};

	inline void assign(const vector2& p1, const vector2& p2, const vector2& p3, const vector2& p4) {
		this->X_i = p1.X; this->Y_i = p1.Y;
		this->X_j = p2.X; this->Y_j = p2.Y;
		this->X_k = p3.X; this->Y_k = p3.Y;
		this->X_l = p4.X; this->Y_l = p4.Y;
	};

public:

	//**********************************************************************************
	vector8() : X_i{}, Y_i{}, X_j{}, Y_j{}, X_k{}, Y_k{}, X_l{}, Y_l{} {}
	vector8(double xi, double yi, double xj, double yj, double xk, double yk, double xl, double yl) :
	X_i{}, Y_i{}, X_j{}, Y_j{}, X_k{}, Y_k{}, X_l{}, Y_l{} { assign(xi,yi,xj,yj,xk,yk,xl,yl); }
	vector8(double a) : X_i{}, Y_i{}, X_j{}, Y_j{}, X_k{}, Y_k{}, X_l{}, Y_l{} { assign(a); };
	vector8(vector8& vec) : X_i{}, Y_i{}, X_j{}, Y_j{}, X_k{}, Y_k{}, X_l{}, Y_l{} { assign(vec); };
	vector8(vector2& p1, vector2& p2, vector2& p3, vector2& p4) :X_i{}, Y_i{}, X_j{}, Y_j{}, X_k{}, Y_k{}, X_l{}, Y_l{}
	{assign(p1, p2, p3, p4);}

	//**********************************************************************************

	inline double VxV(const vector8& a, const  vector8& b) { return (
		a.X_i * b.X_i + a.Y_i * b.Y_i 
		+ a.X_j * b.X_j + a.Y_j * b.Y_j 
		+ a.X_k * b.X_k + a.Y_k * b.Y_k 
		+ a.X_l * b.X_l + a.Y_l * b.Y_l  ); } const
	inline double Norm(const vector8& a) { return std::sqrt(VxV(a, a)); } const
	inline void UniVector() { assign(*this /= Norm(*this)); };
	inline const vector8& UniVector(const vector8& vec) { vector8 res{}; res = vec / Norm(vec); return res; } const

	//**********************************************************************************

	vector8 operator + (const vector8& add) const
	{
		return vector8(
			this->X_i + add.X_i, this->Y_i + add.Y_i,
			this->X_j + add.X_j, this->Y_j + add.Y_j,
			this->X_k + add.X_k, this->Y_k + add.Y_k,
			this->X_l + add.X_l, this->Y_l + add.Y_l);
	}

	vector8 operator - (const vector8& sub) const
	{
		return vector8(
			this->X_i - sub.X_i, this->Y_i - sub.Y_i,
			this->X_j - sub.X_j, this->Y_j - sub.Y_j,
			this->X_k - sub.X_k, this->Y_k - sub.Y_k,
			this->X_l - sub.X_l, this->Y_l - sub.Y_l);
	}

	vector8 operator / (double div) const
	{
		return vector8(
			this->X_i /div, this->Y_i / div,
			this->X_j / div, this->Y_j / div,
			this->X_k / div, this->Y_k / div,
			this->X_l / div, this->Y_l / div);
	}

	vector8 operator * (double mul) const
	{
		return vector8(
			this->X_i * mul, this->Y_i * mul,
			this->X_j * mul, this->Y_j * mul,
			this->X_k * mul, this->Y_k * mul,
			this->X_l * mul, this->Y_l * mul);
	}

	vector8& operator += (const vector8& add)
	{
		assign(*this + add);
		return *this;
	}

	vector8& operator -= (const vector8& sub)
	{
		assign(*this - sub);
		return *this;
	}

	vector8& operator /= (double div)
	{
		assign(*this / div);
		return *this;
	}

	vector8& operator *= (double mul)
	{
		assign(*this * mul);
		return *this;
	}

	vector8& operator = (const vector8& vec)
	{
		assign(vec);
		return *this;
	}
};
//**********************************************************************************
//**********************************************************************************
class vector3 
{

};