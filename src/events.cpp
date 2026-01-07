#include "events.h"

events::events() = default;
events::~events() = default;

void events::add_particle(int pid, double t, double x, double y, double z,
                          double e, double px, double py, double pz)
{
    particles part;
    part.set_pid(pid);
    part.set_px(px);
    part.set_py(py);
    part.set_pz(pz);
    part.set_e(e);
    part.set_t(t);
    part.set_x(x);
    part.set_y(y);
    part.set_z(z);

    particle_vector.push_back(part);
}

void events::add_particle(int pid, double t, double x, double y, double z,
                          double e, double px, double py, double pz,
                          bool recon_flag, double ww)
{
    particles part;
    part.set_pid(pid);
    part.set_px(px);
    part.set_py(py);
    part.set_pz(pz);
    part.set_e(e);
    part.set_t(t);
    part.set_x(x);
    part.set_y(y);
    part.set_z(z);
    part.set_reconst_flag(recon_flag);
    part.set_weight(ww);

    particle_vector.push_back(part);
}

void events::add_particle(const particles& part)
{
    particle_vector.push_back(part);
}


