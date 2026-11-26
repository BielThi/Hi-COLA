// Power spectrum of a Gadget snapshot of a Hi-COLA run, with Hi-COLA's own estimator.
//
//   mpirun -np <n> scripts/pofk_snapshot/pofk_snapshot <gadget file prefix> <boxsize Mpc/h> <npart total>
//          <nmesh> <density assignment method> <interlacing 0|1> <subtract shotnoise 0|1> <output file>
//
// It calls the function Hi-COLA's nbody calls for pofk_z<z>.txt (COLASolver/src/AnalyzeOutput.h,
// compute_power_spectrum -> FML::CORRELATIONFUNCTIONS::compute_power_spectrum), with the same settings:
//   - density assignment NGP | CIC | TSC | PCS | PQS on an nmesh^3 grid, optionally interlaced (a second
//     grid of the particles shifted by half a cell, averaged in Fourier space with the phase e^{i k.shift})
//   - deconvolution of the window, every mode divided by prod_i sinc(pi k_i / 2 k_Ny)^p (p = 2 CIC, 3 TSC,
//     4 PCS), no alias sum
//   - nmesh/2 linear bins of width 2 pi / L from 2 pi / L; k is the mean |k| of the modes of the bin
//   - shot noise 1 / Npart (box units) subtracted after the deconvolution if asked
//   - k in h/Mpc, P(k) in (Mpc/h)^3
// Only the density assignment method changes between calls; with the run's own settings it reproduces
// pofk_z<z>.txt up to the single-precision positions of the snapshot.
//
// Reading: every MPI task reads the files i = task, task + ntasks, ... of the snapshot, and the particles
// are then sent to the task owning their x-slab (the grid's slab decomposition).
#include <FML/ComputePowerSpectra/ComputePowerSpectrum.h>
#include <FML/GadgetUtils/GadgetUtils.h>
#include <FML/Global/Global.h>
#include <FML/MPIParticles/MPIParticles.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

const int NDIM = 3;

// Positions only: the Gadget reader fills the fields the particle has (get_pos), in units of the box [0,1)
struct Particle {
    double pos[NDIM];
    constexpr int get_ndim() const { return NDIM; }
    double * get_pos() { return pos; }
};

using GadgetReader = FML::FILEUTILS::GADGET::GadgetReader;

int main(int argc, char ** argv) {
    if (argc != 9) {
        if (FML::ThisTask == 0)
            std::cerr << "usage: pofk_snapshot <gadget file prefix> <boxsize Mpc/h> <npart total> <nmesh> <method> "
                         "<interlacing 0|1> <subtract shotnoise 0|1> <output file>\n";
        return 1;
    }
    const std::string fileprefix = argv[1];
    const double boxsize = std::stod(argv[2]);
    const size_t npart_expected = std::stoull(argv[3]);
    const int nmesh = std::stoi(argv[4]);
    const std::string method = argv[5];
    const bool interlacing = std::stoi(argv[6]) != 0;
    const bool subtract_shotnoise = std::stoi(argv[7]) != 0;
    const std::string output = argv[8];

    //=============================================================
    // Read the snapshot: this task's share of the files
    //=============================================================
    GadgetReader reader;
    const int nfiles = reader.get_num_files(fileprefix + ".0");
    const double header_boxsize = reader.get_header().BoxSize;
    std::vector<Particle> part;
    for (int i = FML::ThisTask; i < nfiles; i += FML::NTasks) {
        const bool only_keep_part_in_domain = false;
        const bool verbose = false;
        reader.read_gadget_single(fileprefix + "." + std::to_string(i), part, only_keep_part_in_domain, verbose);
    }
    // Positions exactly at 1 (single precision rounding) belong to the periodic image at 0
    for (auto & p : part)
        for (int idim = 0; idim < NDIM; idim++) {
            if (p.pos[idim] >= 1.0)
                p.pos[idim] -= 1.0;
            if (p.pos[idim] < 0.0)
                p.pos[idim] += 1.0;
        }

    // Send every particle to the task of its x-slab; room for 3x the mean number per task
    FML::PARTICLE::MPIParticles<Particle> particles;
    const size_t nmean = npart_expected / FML::NTasks + 1;
    FML::Vector<Particle> local;
    local.reserve(std::max(3 * nmean, part.size()));
    local.assign(part.begin(), part.end());
    std::vector<Particle>().swap(part);
    particles.move_from(std::move(local));
    particles.communicate_particles();

    const size_t npart = particles.get_npart_total();
    if (FML::ThisTask == 0) {
        std::cout << "Snapshot:        " << fileprefix << ".* (" << nfiles << " files, header BoxSize " << header_boxsize
                  << ")\n";
        std::cout << "Particles:       " << npart << " (expected " << npart_expected << ")\n";
        std::cout << "Box:             " << boxsize << " Mpc/h\n";
        std::cout << "P(k):            nmesh " << nmesh << ", " << method << ", interlacing " << interlacing
                  << ", shot noise subtracted " << subtract_shotnoise << "\n";
        std::cout << "MPI tasks:       " << FML::NTasks << "\n";
    }
    if (npart != npart_expected) {
        if (FML::ThisTask == 0)
            std::cerr << "Read " << npart << " particles, expected " << npart_expected << "\n";
        return 1;
    }

    //=============================================================
    // P(k) exactly as in Hi-COLA's compute_power_spectrum
    //=============================================================
    FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<NDIM> pofk(nmesh / 2);
    pofk.subtract_shotnoise = subtract_shotnoise;
    FML::CORRELATIONFUNCTIONS::compute_power_spectrum<NDIM, Particle>(
        nmesh, particles.get_particles_ptr(), particles.get_npart(), npart, pofk, method, interlacing);
    pofk.scale(boxsize);

    if (FML::ThisTask == 0) {
        std::ofstream fp(output);
        fp << "# P(k) of " << fileprefix << ".* with Hi-COLA's estimator (scripts/pofk_snapshot)\n";
        fp << "# nmesh = " << nmesh << ", density assignment = " << method << ", interlacing = " << std::boolalpha
           << interlacing << ", ShotnoiseSubtracted = " << subtract_shotnoise << ", npart = " << npart
           << ", boxsize = " << boxsize << " Mpc/h\n";
        fp << "# k  (h/Mpc)    P(k)  (Mpc/h)^3\n";
        for (int i = 0; i < pofk.n; i++)
            fp << std::setw(15) << pofk.kbin[i] << " " << std::setw(15) << pofk.pofk[i] << "\n";
        std::cout << "Wrote " << output << "\n";
    }
    return 0;
}
