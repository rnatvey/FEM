#pragma once
class ElementData
{

private:

	double E{};
	double mu{};
	double G{ E / (1.0 - 2.0 * mu) };
	size_t id{ 0 };

public:

	ElementData() : E{}, mu{} {}
	ElementData(double E, double mu) : E{ E }, mu{ mu } {}

	//auto get_C() const -> decltype(std::decval<Derived>().C);

	inline double get_E() const { return E; } 
	inline double get_mu() const { return mu; } 
	inline double get_G() const { return G; } 
	inline size_t get_id() const { return id; } 

	inline void assign_E(double E) { this->E = E; }
	inline void assign_mu(double mu) { this->mu = mu; }
	inline void assign_G(double G) { this->G = G; }
	inline void assign_id(size_t id) { this->id = id; }

};