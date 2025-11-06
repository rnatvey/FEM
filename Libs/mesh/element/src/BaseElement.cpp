#include "BaseElement.h"
#include "node.h"
#include "material.h"

BaseElement::BaseElement(int id, const std::vector<int>& nodeIds, int materialId)
    : id_(id), nodeIds_(nodeIds), materialId_(materialId) {
}