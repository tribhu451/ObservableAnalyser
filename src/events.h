#pragma once
#include <vector>
#include "particles.h"

class events {
public:
    events();
    ~events();

    void add_particle(int pid, double t, double x, double y,
                      double z, double e, double px, double py, double pz);

    void add_particle(int pid, double t, double x, double y,
                      double z, double e, double px, double py, double pz,
                      bool recon_flag, double weight);

    void add_particle(const particles& part);

    particles* get_particle(int i) {
        return &particle_vector[i];   // OK: pointer to internal object
    }

    int get_multiplicity_of_the_event() const {
        return particle_vector.size();
    }

    void clear_particle_vector() {
        particle_vector.clear();
    }

private:
    std::vector<particles> particle_vector;
};


