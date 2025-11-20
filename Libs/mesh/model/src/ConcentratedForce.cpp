#include "ConcentratedForce.h"

ConcentratedForce::ConcentratedForce(int nodeId, double fx, double fy)
    : nodeId_(nodeId), fx_(fx), fy_(fy), force_(fx, fy) {
}