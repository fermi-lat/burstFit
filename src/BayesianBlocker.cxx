/** \file BayesianBlocker.cxx
    \brief Implementation of BayesianBlocker class.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "CLHEP/Random/Stat.h"

#include "burstFit/BayesianBlocker.h"

namespace burstFit {

  const int BayesianBlocker::s_ncp_prior;

  BayesianBlocker::BayesianBlocker(double cell_size, const std::vector<double> & cell_pop): m_cell_size(cell_pop.size(), cell_size),
    m_cell_pop(cell_pop), m_blocks(), m_blocks_computed(false) {
  }

  const std::vector<double> & BayesianBlocker::getBlocks() {
    if (!m_blocks_computed) {
      computeBlocks();
      m_blocks_computed = true;
    }
    return m_blocks;
  }

  void BayesianBlocker::computeBlocks() {
    // Number of cells, obtained here once for convenience.
    vec_t::size_type num_cells = m_cell_pop.size();

    // Create arrays to hold the reverse cumulative sums.
    deque_t rev_csize(1, m_cell_size[0]);
    deque_t rev_cpop(1, m_cell_pop[0]);

    vec_t best(num_cells, 0.);
    vec_t merged(num_cells, 0.);
    vec_t temp(num_cells, 0.);
    std::vector<vec_t::size_type> last_start(num_cells, 0);
    for (vec_t::size_type index = 0; index != num_cells; ++index) {

      if (index > 0) {
        // The following lines replace the following construct from BBglobal.pro:
        //       cumsizes = [cell_sizes(R), cumsizes + cell_sizes(R)]
        for (deque_t::iterator itor = rev_csize.begin(); itor != rev_csize.begin() + index; ++itor) *itor += m_cell_size[index];
        rev_csize.push_front(m_cell_size[index]);

        // The following lines replace the following construct from BBglobal.pro:
        //       cumpops  = [cell_pops(R),  cumpops  + cell_pops(R)]
        for (deque_t::iterator itor = rev_cpop.begin(); itor != rev_cpop.begin() + index; ++itor) *itor += m_cell_pop[index];
        rev_cpop.push_front(m_cell_pop[index]);
      }

      // The following lines replace the following construct from BBglobal.pro:
      //    merged = reverse(log_prob(cumpops, cumsizes))
      computeLogProb(rev_csize, rev_cpop, merged);
      std::reverse(merged.begin(), merged.begin() + index + 1);

      // The following lines replace the following construct from BBglobal.pro:
      //    if (R eq 0) then begin
      //       best(R) = max(merged, imaxer)
      //     endif else begin
      //       temp = [0., best(0:R-1)] + merged
      //       best(R) = max(temp, imaxer)
      //     endelse
      vec_t::size_type imaxer;
      vec_t::const_iterator imax_itor;
      if (index == 0) {
        imax_itor = std::max_element(merged.begin(), merged.begin() + index + 1);
        imaxer = imax_itor - merged.begin();
      } else {
        temp[0] = 0.;
        std::copy(best.begin(), best.begin() + index, temp.begin() + 1);
        for (vec_t::size_type ii = 0; ii != index + 1; ++ii) temp[ii] += merged[ii];
        imax_itor = std::max_element(temp.begin(), temp.begin() + index + 1);
        imaxer = imax_itor - temp.begin();
      }
      best[index] = *imax_itor;

      // The following lines replace the following construct from BBglobal.pro:
      //    last_start(R) = imaxer
      last_start[index] = imaxer;
    }

    // The following lines replace the following construct from BBglobal.pro:
    // ; Find the optimum partition, thereby generating the changepoint indices (CPs).
    //
    // CPs = long(num_cells)                           ;initialize CPs array with end point of time series
    // index = last_start(num_cells-1)                 ;the last-found CP
    //
    // while (index gt 1) do begin
    //    CPs = [index, CPs]                           ;fish out from last_start the CPs found above
    //    index = last_start(index-1)
    //  endwhile
    //
    // CPs = [0, CPs]                                  ;finalize the CPs array with start point of time series
    std::deque<vec_t::size_type> cp(1, num_cells);
    for (vec_t::size_type index = last_start[num_cells - 1]; index > 1; index = last_start[index - 1])
      cp.push_front(index);
    cp.push_front(0);

  std::cout << "CPs: " << std::endl;
  for (std::deque<vec_t::size_type>::iterator itor = cp.begin(); itor != cp.end(); ++itor) std::cout << "\t" << *itor;
  std::cout << std::endl;
  }

  void BayesianBlocker::computeLogProb(const std::deque<double> & rev_csize, const std::deque<double> & rev_cpop,
    std::vector<double> & log_prob) const {
    // This code is for data type "3" from pulsefitter6.pro.
    //
    // The following lines replace the following construct from log_prob.pro:
    // ncells = n_elements(cell_sizes)
    deque_t::size_type num_cells = rev_csize.size();

    // The following lines replace the following construct from log_prob.pro:
    // ;-----------------------------
    // ; posterior for binned data
    // ;-----------------------------
    //
    //    logprob = lngamma(cell_pops + 1.) - (cell_pops + 1.) * alog(cell_sizes)
    //
    //    logprob = logprob - ncp_prior
    for (deque_t::size_type index = 0; index != num_cells; ++index)
      log_prob[index] = HepStat::gammln(rev_cpop[index] + 1.) - (rev_cpop[index] + 1.) * log(rev_csize[index]) - s_ncp_prior;
  }

}
