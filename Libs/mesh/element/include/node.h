#pragma once
#include <Eigen/Dense>

class Node {
public:
    Node(int id, double x, double y) : id_(id), coordinates_(x, y) {}

    int getId() const { return id_; }
    const Eigen::Vector2d& getCoordinates() const { return coordinates_; }
    void setCoordinates(const Eigen::Vector2d& coords) { coordinates_ = coords; }

private:
    int id_;
    Eigen::Vector2d coordinates_;
};