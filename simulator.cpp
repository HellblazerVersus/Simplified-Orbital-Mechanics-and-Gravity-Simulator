#include <iostream>
#include <fstream>
#include <string>

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

template <typename T>
T maxVal(T a, T b) { return a > b ? a : b; }

class CollisionError {
    std::string message;
public:
    explicit CollisionError(const std::string& msg) : message(msg) {}
    const char* what() const noexcept { return message.c_str(); } 
};

class SpaceObject {
public:
    static constexpr double G = 6.67430e-11; 

    virtual void displayStatus() const = 0; 
    virtual ~SpaceObject() {}
};


class CelestialBody;
void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time);

class CelestialBody : public SpaceObject {
private:
    double mass;
    double x, y;
    double vx, vy;

public:
    inline double getX() const { return x; }
    inline double getY() const { return y; }
    inline double getMass() const { return mass; }


    CelestialBody()
        : mass(0.0), x(0.0), y(0.0), vx(0.0), vy(0.0) {}

    
    CelestialBody(double m, double px, double py, double pvx = 0.0, double pvy = 0.0)
        : mass(m), x(px), y(py), vx(pvx), vy(pvy) {}


    ~CelestialBody() {
      
    }

    void calculateForce(const CelestialBody& other, double& fx, double& fy) const {
        double dx = other.x - this->x;
        double dy = other.y - this->y;
        double r2 = sqr(dx) + sqr(dy);
        if (r2 <= 0.0) {
            
            throw CollisionError("Collision or zero distance detected between bodies.");
        }
        double r = sqrt_newton(r2);
        double F = (SpaceObject::G * this->mass * other.mass) / r2;
        fx = F * (dx / r);
        fy = F * (dy / r);
    }

    void updateVelocity(double fx, double fy, double dt) {
        
        this->vx += (fx / this->mass) * dt;
        this->vy += (fy / this->mass) * dt;
    }

    void updatePosition(double dt) {
        x += vx * dt;
        y += vy * dt;
    }

    void displayStatus() const override {
        std::cout << "[CelestialBody] mass=" << mass << " pos=(" << x << "," << y << ") vel=(" << vx << "," << vy << ")\n";
    }


    friend void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time);


    friend std::ostream& operator<<(std::ostream& os, const CelestialBody& b) {
        os << "(" << b.x << ", " << b.y << ")";
        return os;
    }
};


void logPosition(const CelestialBody& body, std::ofstream& ofs, const std::string& id, double time) {
    ofs << id << "," << time << "," << body.x << "," << body.y << "\n";
}


void display(const CelestialBody& b) {
    std::cout << "Display body at " << b << "\n";
}
void display(const std::string& s) {
    std::cout << "Display message: " << s << "\n";
}


void reportStatus(const SpaceObject* obj) {
    obj->displayStatus();
}

int main() {

    const double AU = 1.496e11; 
    const double sunMass = 1.98847e30;
    const double earthMass = 5.972e24;

    CelestialBody bodies[2] = {
        CelestialBody(sunMass, 0.0, 0.0, 0.0, 0.0),                   
        CelestialBody(earthMass, AU, 0.0, 0.0, 29780.0)                
    };

  
    CelestialBody* dynamicPlanet = new CelestialBody(earthMass, AU * 1.05, 0.0, 0.0, 29000.0);


    double largerMass = maxVal(bodies[0].getMass(), dynamicPlanet->getMass());
    std::cout << "Larger mass between Sun and dynamic planet: " << largerMass << "\n";

 
    std::ofstream traj("trajectory.csv");
    if (!traj) {
        std::cerr << "Failed to open trajectory.csv for writing\n";
        delete dynamicPlanet;
        return 1;
    }
    traj << "id,time,x,y\n";

    
    reportStatus(&bodies[0]);          
    reportStatus(&bodies[1]);          
    reportStatus(dynamicPlanet);     

    
    double dt = 60.0; 
    int steps = 24 * 60; 

    
    double simTime = 0.0;
    logPosition(bodies[1], traj, "planet", simTime);
    logPosition(*dynamicPlanet, traj, "dynamic", simTime);


    for (int step = 0; step < steps; ++step) {
        simTime += dt;
       
        try {
            double fx12 = 0.0, fy12 = 0.0;
            bodies[1].calculateForce(bodies[0], fx12, fy12);
         
            bodies[1].updateVelocity(fx12, fy12, dt);
            bodies[1].updatePosition(dt);

            double fx_dyn = 0.0, fy_dyn = 0.0;
            dynamicPlanet->calculateForce(bodies[0], fx_dyn, fy_dyn);
            dynamicPlanet->updateVelocity(fx_dyn, fy_dyn, dt);
            dynamicPlanet->updatePosition(dt);

       
            double fx_sun = -fx12 - fx_dyn;
            double fy_sun = -fy12 - fy_dyn;
            bodies[0].updateVelocity(fx_sun, fy_sun, dt);
            bodies[0].updatePosition(dt);
        }
        catch (const CollisionError& e) {
            std::cerr << "Simulation error: " << e.what() << " at step " << step << "\n";
            break;
        }

        logPosition(bodies[1], traj, "planet", simTime);
        logPosition(*dynamicPlanet, traj, "dynamic", simTime);
    }


    display(bodies[1]);
    display(std::string("End of simulation"));


    std::cout << "Final planet position: " << bodies[1] << "\n";
    std::cout << "Final dynamic planet position: " << *dynamicPlanet << "\n";


    delete dynamicPlanet;

    traj.close();

    return 0;

}
