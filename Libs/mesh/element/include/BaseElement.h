#pragma once
#include "matrix.h"
#include "ElementData.h"

class BaseElement {

private:

	ElementData& data;

public:

	virtual auto get_K() const = 0; //Возвращает матрицу жесткости в глобальной системе координат
	virtual auto get_f() const = 0; //Возвращает перемещения в глобальной системе координат

};
