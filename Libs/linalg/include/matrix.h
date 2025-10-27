#pragma once
#include "vector.h"

class matrix2 
{

};
//**********************************************************************************
//**********************************************************************************
class matrix3 
{
};
//**********************************************************************************
//**********************************************************************************
class matrix3x8 
{

};
//**********************************************************************************
//**********************************************************************************
class matrix8 
{
public:
	
	enum AxisIndex
	{
		Xi,Yi,
		Xj,Yj,
		Xk,Yk,
		Xl,Yl
	};

	double xixi{ 0.0 }, xiyi{ 0.0 }, xixj{ 0.0 }, xiyj{ 0.0 }, xixk{ 0.0 }, xiyk{ 0.0 }, xixl{ 0.0 }, xiyl{ 0.0 };
	double yixi{ 0.0 }, yiyi{ 0.0 }, yixj{ 0.0 }, yiyj{ 0.0 }, yixk{ 0.0 }, yiyk{ 0.0 }, yixl{ 0.0 }, yiyl{ 0.0 };
	double xjxi{ 0.0 }, xjyi{ 0.0 }, xjxj{ 0.0 }, xjyj{ 0.0 }, xjxk{ 0.0 }, xjyk{ 0.0 }, xjxl{ 0.0 }, xjyl{ 0.0 };
	double yjxi{ 0.0 }, yjyi{ 0.0 }, yjxj{ 0.0 }, yjyj{ 0.0 }, yjxk{ 0.0 }, yjyk{ 0.0 }, yjxl{ 0.0 }, yjyl{ 0.0 };
	double xkxi{ 0.0 }, xkyi{ 0.0 }, xkxj{ 0.0 }, xkyj{ 0.0 }, xkxk{ 0.0 }, xkyk{ 0.0 }, xkxl{ 0.0 }, xkyl{ 0.0 };
	double ykxi{ 0.0 }, ykyi{ 0.0 }, ykxj{ 0.0 }, ykyj{ 0.0 }, ykxk{ 0.0 }, ykyk{ 0.0 }, ykxl{ 0.0 }, ykyl{ 0.0 };
	double xlxi{ 0.0 }, xlyi{ 0.0 }, xlxj{ 0.0 }, xlyj{ 0.0 }, xlxk{ 0.0 }, xlyk{ 0.0 }, xlxl{ 0.0 }, xlyl{ 0.0 };
	double ylxi{ 0.0 }, ylyi{ 0.0 }, ylxj{ 0.0 }, ylyj{ 0.0 }, ylxk{ 0.0 }, ylyk{ 0.0 }, ylxl{ 0.0 }, ylyl{ 0.0 };

private:
	void assignByRow(int axInd, const vector8& ax);
	void assignByCol(int axInd, const vector8& ax);
	inline void assign(const matrix8& matrix) { return; }
	inline void assign(double a) { return; }
	inline void assign(int i, int j, double elem) { return; }
	inline void MxV(const matrix8& mat, const vector8& vec, vector8& result) { return; }
	inline void MxM(const matrix8& mat1, const matrix8& mat2, matrix8& result) { return; }
};