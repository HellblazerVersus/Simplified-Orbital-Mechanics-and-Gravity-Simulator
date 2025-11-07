#include <iostream>
#include <fstream>
#include <string>

// Simple math helpers (no <cmath> allowed)
inline double absd(double v) { return v < 0.0 ? -v : v; }
double sqrt_newton(double n) {
    if (n <= 0.0) return 0.0;
    double x = n;
    for (int i = 0; i < 40; ++i) {
        double nx = 0.5 * (x + n / x);
        if (absd(nx - x) < 1e-12) break;
        x = nx;
    }
    return x;
}
inline double sqr(double v) { return v * v; }

// Template utility (SN 5)
template <typename T>
T maxVal(T a, T b) { return a > b ? a : b; }

// Custom exception for collisions / division by zero scenario (SN 5)
class CollisionError {
    std::string message;
public:
    explicit CollisionError(const std::string& msg) : message(msg) {}
    const char* what() const noexcept { return message.c_str(); } // made noexcept
};

// Abstract base class (SN 3)
class SpaceObject {
public:
    static constexpr double G = 6.67430e-11; // static constant gravitational constant (SN 4)

    virtual void displayStatus() const = 0; // pure virtual -> abstract
    virtual ~SpaceObject() {}
};

// Forward declaration for friend function
class CelestialBody;
void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time);

// Derived class (SN 1, SN 3, SN 4)
class CelestialBody : public SpaceObject {
private:
    double mass;
    double x, y;
    double vx, vy;

public:
    // Inline getters (SN 2)
    inline double getX() const { return x; }
    inline double getY() const { return y; }
    inline double getMass() const { return mass; }

    // Default constructor
    CelestialBody()
        : mass(0.0), x(0.0), y(0.0), vx(0.0), vy(0.0) {}

    // Custom constructor (SN 2)
    CelestialBody(double m, double px, double py, double pvx = 0.0, double pvy = 0.0)
        : mass(m), x(px), y(py), vx(pvx), vy(pvy) {}

    // Destructor (SN 2)
    ~CelestialBody() {
        // minimal cleanup (no dynamic members here), just illustrative
    }

    // calculateForce: computes gravitational force exerted by 'other' on this (SN 1)
    // returns components via reference parameters
    void calculateForce(const CelestialBody& other, double& fx, double& fy) const {
        double dx = other.x - this->x;
        double dy = other.y - this->y;
        double r2 = sqr(dx) + sqr(dy);
        if (r2 <= 0.0) {
            // potential division by zero: throw custom exception (SN 5)
            throw CollisionError("Collision or zero distance detected between bodies.");
        }
        double r = sqrt_newton(r2);
        double F = (SpaceObject::G * this->mass * other.mass) / r2;
        fx = F * (dx / r);
        fy = F * (dy / r);
    }

    // updateVelocity uses 'this' pointer explicitly (SN 2)
    void updateVelocity(double fx, double fy, double dt) {
        // this pointer used explicitly
        this->vx += (fx / this->mass) * dt;
        this->vy += (fy / this->mass) * dt;
    }

    // Update position using current velocity
    void updatePosition(double dt) {
        x += vx * dt;
        y += vy * dt;
    }

    // Overridden displayStatus (SN 3)
    void displayStatus() const override {
        std::cout << "[CelestialBody] mass=" << mass << " pos=(" << x << "," << y << ") vel=(" << vx << "," << vy << ")\n";
    }

    // Friend function to allow logging access to private members (SN 2)
    friend void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time);

    // Stream insertion operator overload (SN 4)
    friend std::ostream& operator<<(std::ostream& os, const CelestialBody& b) {
        os << "(" << b.x << ", " << b.y << ")";
        return os;
    }
};

// Friend function definition (SN 2) - logs pos to CSV
void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time) {
    ofs << id << "," << time << "," << body.x << "," << body.y << "\n";
}

// Function overloading for display (SN 2)
void display(const CelestialBody& b) {
    std::cout << "Display body at " << b << "\n";
}
void display(const std::string& s) {
    std::cout << "Display message: " << s << "\n";
}

// Demonstrate polymorphism and dynamic binding (SN 4)
void reportStatus(const SpaceObject* obj) {
    obj->displayStatus();
}

int main() {
    // Create an array of CelestialBody objects size 2 (SN 1)
    // Sun at origin, Planet at 1 AU on x-axis with a tangential velocity for approximate circular orbit
    const double AU = 1.496e11; // meters
    const double sunMass = 1.98847e30;
    const double earthMass = 5.972e24;

    CelestialBody bodies[2] = {
        CelestialBody(sunMass, 0.0, 0.0, 0.0, 0.0),                      // Sun
        CelestialBody(earthMass, AU, 0.0, 0.0, 29780.0)                 // Planet
    };

    // Dynamic allocation of a CelestialBody (SN 2)
    CelestialBody* dynamicPlanet = new CelestialBody(earthMass, AU * 1.05, 0.0, 0.0, 29000.0);

    // Demonstrate template usage
    double largerMass = maxVal(bodies[0].getMass(), dynamicPlanet->getMass());
    std::cout << "Larger mass between Sun and dynamic planet: " << largerMass << "\n";

    // Open file for trajectory logging (SN 5)
    std::ofstream traj("trajectory.csv");
    if (!traj) {
        std::cerr << "Failed to open trajectory.csv for writing\n";
        delete dynamicPlanet;
        return 1;
    }
    traj << "id,time,x,y\n";

    // Polymorphism: use base-class pointers (SN 4)
    reportStatus(&bodies[0]);           // Sun
    reportStatus(&bodies[1]);           // Planet
    reportStatus(dynamicPlanet);        // dynamic planet

    // Simulation parameters
    double dt = 60.0; // seconds per step
    int steps = 24 * 60; // simulate one day in 1-minute steps

    // Write initial positions (time = 0)
    double simTime = 0.0;
    logPosition(bodies[1], traj, "planet", simTime);
    logPosition(*dynamicPlanet, traj, "dynamic", simTime);

    // Simulation loop (very simple Euler integrator)
    for (int step = 0; step < steps; ++step) {
        simTime += dt;
        // Compute force between sun and planet (bodies[1] due to bodies[0])
        try {
            double fx12 = 0.0, fy12 = 0.0;
            bodies[1].calculateForce(bodies[0], fx12, fy12);
            // Update planet velocity and position
            bodies[1].updateVelocity(fx12, fy12, dt);
            bodies[1].updatePosition(dt);

            // Also update dynamicPlanet affected by sun
            double fx_dyn = 0.0, fy_dyn = 0.0;
            dynamicPlanet->calculateForce(bodies[0], fx_dyn, fy_dyn);
            dynamicPlanet->updateVelocity(fx_dyn, fy_dyn, dt);
            dynamicPlanet->updatePosition(dt);

            // Basic mutual reaction: accelerate sun a tiny bit from planets (Newton's 3rd law)
            double fx_sun = -fx12 - fx_dyn;
            double fy_sun = -fy12 - fy_dyn;
            bodies[0].updateVelocity(fx_sun, fy_sun, dt);
            bodies[0].updatePosition(dt);
        }
        catch (const CollisionError& e) {
            std::cerr << "Simulation error: " << e.what() << " at step " << step << "\n";
            break;
        }

        // Log positions each step (SN 5 via friend function)
        logPosition(bodies[1], traj, "planet", simTime);
        logPosition(*dynamicPlanet, traj, "dynamic", simTime);
    }

    // Demonstrate overloaded display functions (SN 2)
    display(bodies[1]);
    display(std::string("End of simulation"));

    // Use overloaded operator<< (SN 4)
    std::cout << "Final planet position: " << bodies[1] << "\n";
    std::cout << "Final dynamic planet position: " << *dynamicPlanet << "\n";

    // Clean up dynamic allocation (SN 2)
    delete dynamicPlanet;

    traj.close();

    return 0;
}