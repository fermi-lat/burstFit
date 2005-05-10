/** \file BayesianBinner.h
    \brief Declaration of BayesianBinner class.
*/
#ifndef evtbin_BayesianBinner_h
#define evtbin_BayesianBinner_h

#include <deque>
#include <string>
#include <vector>

#include "evtbin/Binner.h"

namespace evtbin {
  /** \class BayesianBinner
      \brief Declaration of a linearly uniform interval binner.
  */
  class BayesianBinner : public Binner {
    public:
      typedef std::deque<double> deque_t;
      typedef std::vector<double> vec_t;

      /** \brief Construct a linear binner object.
          \param interval_begin Left boundary of the binning interval.
          \param interval_end Right boundary of the binning interval.
          \param bin_size The size of the bins.
          \param name Optional name of the quantity being binned.
      */
      template <typename Itor>
      BayesianBinner(const std::string & name, double cell_size, Itor begin, Itor end): Binner(name),
        m_cell_size(end - begin, cell_size), m_cell_pop(begin, end), m_cells() { computeBlocks(); }

      virtual ~BayesianBinner() throw();

      /** \brief Return the bin number for the given value.
          \param value The value being binned.
      */
      virtual long computeIndex(double value) const;

      /** \brief Return the number of bins currently defined.
      */
      virtual long getNumBins() const;

      /** \brief Return the interval spanned by the given bin.
          \param index The index indicating the bin number.
      */
      virtual Binner::Interval getInterval(long index) const;

      /** \brief Create copy of this object.
      */
      virtual Binner * clone() const;

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
      std::vector<double> m_cells;
  };

}

#endif
