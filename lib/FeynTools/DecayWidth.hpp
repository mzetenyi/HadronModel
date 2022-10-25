#ifndef DECAYWIDTH_HPP
#define DECAYWIDTH_HPP

#include "utils.hpp"

/**
 * @brief Represents the decay width of a resonance decaying into two stable
 * particles.
 *
 *
 * @tparam Amplitude An object representing the amplitude of the decay process.
 * Should have a constructor DecayWidthSS(double M) taking the resonance mass as
 * amplitude, and a method double MSQR() returning the squared amplitude
 * (summed/averaged over polarizations).
 */
template <typename Amplitude>
class DecayWidthSS {
 public:
  DecayWidthSS(Amplitude amp) : amp(amp) {}
  double GammaTot() {
    double p = amp.getKinematics().getpabs();
    double M = amp.M;
    return 1. / (8. * pi_) * amp.MSQR() * p / (M * M);
  }

 private:
  Amplitude amp;
};

#endif  // DECAYWIDTH_HPP