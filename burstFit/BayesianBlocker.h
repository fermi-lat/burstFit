/** \file BayesianBlocker.h
    \brief Declaration for BayesianBlocker class.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#ifndef burstFit_BayesianBlocker_h
#define burstFit_BayesianBlocker_h

#include <deque>
#include <vector>

namespace burstFit {

  class BayesianBlocker {
    public:
      typedef std::deque<double> deque_t;
      typedef std::vector<double> vec_t;

      /** \brief Create blocker with constant cell sizes for the given population.
          \param cell_size The (constant) size for each cell.
          \param cell_pop The population (counts) for the data series.
      */
      BayesianBlocker(double cell_size, const std::vector<double> & cell_pop);

      /** \brief Compute Bayesian Blocks if needed, and return them.
      */
      const std::vector<double> & getBlocks();

    private:
      static const int s_ncp_prior = 6;
      /** \brief Perform the Bayesian Block procedure to determine the block definitions.
      */
      void computeBlocks();

      /** \brief Internal utility to compute the log posterior (Bayes factor).
      */
      void computeLogProb(const std::deque<double> & rev_csize, const std::deque<double> & rev_cpop,
        std::vector<double> & result) const;

      std::vector<double> m_cell_size;
      std::vector<double> m_cell_pop;
      std::vector<double> m_blocks;
      bool m_blocks_computed;
  };

}

#endif
