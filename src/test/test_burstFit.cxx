/** \file test_burstFit.cxx
    \brief Test application for burstFit.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#include <iostream>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Random/Stat.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Env.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

class TestBurstFitApp : public st_app::StApp {
  public:
    virtual ~TestBurstFitApp() throw() {}

    virtual void run();

    virtual void testFile(const std::string & file_name);

  private:
    std::string m_data_dir;
};

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

const int BayesianBlocker::s_ncp_prior;

BayesianBlocker::BayesianBlocker(double cell_size, const std::vector<double> & cell_pop): m_cell_size(cell_pop.size(), cell_size),
  m_cell_pop(cell_pop), m_blocks(), m_blocks_computed(false) {
////////////////////////////////////////////////////////////////////////////////
// From make_cells.pro:
////////////////////////////////////////////////////////////////////////////////
//
// if (datatype eq 3) then begin                       ;the Poisson-distributed data are binned:  counts in equal time intervals
//    cell_data = dblarr(num_cells, 2)
//    cellsize = microfact*binscale
//    cell_data(*,0) = cellsize
//    cell_data(*,1) = dat
//  endif
////////////////////////////////////////////////////////////////////////////////
}

const std::vector<double> & BayesianBlocker::getBlocks() {
  if (!m_blocks_computed) {
    computeBlocks();
    m_blocks_computed = true;
  }
  return m_blocks;
}

void BayesianBlocker::computeBlocks() {
////////////////////////////////////////////////////////////////////////////////
// From BBglobal.pro:
////////////////////////////////////////////////////////////////////////////////
//
// time1 = systime(1)
// 
// ; extract the volume and count arrays from cell_data
// 
// cell_sizes = cell_data(*,0)
// cell_pops  = cell_data(*,1)
// 
// num_cells = n_elements(cell_sizes)
// 
// best = dblarr(num_cells)                        ;best(R) will be the value of the optimum at iteration R.
// last_start = lonarr(num_cells)                  ;last(R) will be the index at which this optimum occurs.
// 
// cumsizes = cell_sizes(0)
// cumpops  = cell_pops(0)
// 
// print, f="('ncells = ',I7)", num_cells
// print, 'R = '
// 
// for R = 0L,num_cells-1L do begin
// 
//    if ((R MOD 500) eq 0) then print, f="(I6,$)", R
//    if ((R MOD 10000) eq 0) then print, ''
// 
//    if (R gt 0) then begin
//       cumsizes = [cell_sizes(R), cumsizes + cell_sizes(R)]
//       cumpops  = [cell_pops(R),  cumpops  + cell_pops(R)]
//     endif
// 
//    merged = reverse(log_prob(cumpops, cumsizes))
// 
//    if (R eq 0) then begin
//       best(R) = max(merged, imaxer)
//     endif else begin
//       temp = [0., best(0:R-1)] + merged
//       best(R) = max(temp, imaxer)
//     endelse
// 
//    last_start(R) = imaxer
// 
//  endfor
// 
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
// 
// print, ' '
// print, 'CPs:'
// print, f="(15I8)", CPs
// 
// time2 = systime(1)
// timeCellopt = timeCellopt + time2 - time1
// 
// savbest = best
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// From blocker.pro:
////////////////////////////////////////////////////////////////////////////////
//
// ; Compute block boundaries (tt_blocks) & counts/{pltscale bin} (yy_blocks)
// 
// if (datatype lt 3) then begin                ;time-tagged event data
// 
//    tt_blocks = data(CPs) / microfact         ;convert the CPs to seconds
//    num_blocks = n_elements(tt_blocks) - 1
//    nblocks(itrig) = num_blocks
// 
//    dts = tt_blocks(1:num_blocks) - tt_blocks(0:num_blocks-1)
//    ncounts = nspill * (CPs(1:num_blocks) - CPs(0:num_blocks-1))
//    yy_blocks = pltscale * ncounts / dts      ;intensities in counts/plotting bin
//    yy_blocks = [yy_blocks, 0.]
// 
//  endif else begin                            ;binned data
// 
//    tt_blocks = CPs*binscale                  ;convert the CPs to seconds
//    num_blocks = n_elements(tt_blocks) - 1
//    nblocks(itrig) = num_blocks
//    yy_blocks = fltarr(num_blocks)
//    pltfact = pltscale/binscale
// 
//    for iblock = 0,num_blocks-1 do begin      ;intensities in counts/plotting bin
//       yy_blocks(iblock) = total( data(CPs(iblock):CPs(iblock+1)-1) )
//       yy_blocks(iblock) = pltfact * yy_blocks(iblock) / (CPs(iblock+1) - CPs(iblock))
//     endfor
//    yy_blocks = [yy_blocks, 0.]
// 
//  endelse
////////////////////////////////////////////////////////////////////////////////

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
      // The following lines replace these:
      //       cumsizes = [cell_sizes(R), cumsizes + cell_sizes(R)]
      for (deque_t::iterator itor = rev_csize.begin(); itor != rev_csize.begin() + index; ++itor) *itor += m_cell_size[index];
      rev_csize.push_front(m_cell_size[index]);

      //       cumpops  = [cell_pops(R),  cumpops  + cell_pops(R)]
      for (deque_t::iterator itor = rev_cpop.begin(); itor != rev_cpop.begin() + index; ++itor) *itor += m_cell_pop[index];
      rev_cpop.push_front(m_cell_pop[index]);
    }

    //    merged = reverse(log_prob(cumpops, cumsizes))
    computeLogProb(rev_csize, rev_cpop, merged);
    std::reverse(merged.begin(), merged.begin() + index + 1);

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

//    last_start(R) = imaxer
    last_start[index] = imaxer;
  }

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
////////////////////////////////////////////////////////////////////////////////
// From log_prob.pro
////////////////////////////////////////////////////////////////////////////////
// eps = 1.E-38
  const double eps = std::numeric_limits<float>::epsilon();

// ncells = n_elements(cell_sizes)
  deque_t::size_type num_cells = rev_csize.size();

// if (datatype lt 3) then begin                   ;datatype = 1 or datatype = 2
// 
// ;----------------------------------------
// ; posterior for time-tagged event data
// ;----------------------------------------
// 
//    arg = cell_sizes - cell_pops + 1
//    ii = where (arg gt 0, ngt0)
//    num_bad = ncells - ngt0
// 
//    if (num_bad eq 0) then begin
//       logprob = lngamma(cell_pops + 1.) + lngamma(arg) - lngamma(cell_sizes + 2.)
//     endif else begin
//       logprob = eps + fltarr(ncells)
//       logprob(ii) = lngamma(cell_pops(ii) + 1.) + lngamma(arg(ii)) - lngamma(cell_sizes(ii) + 2.)
//     endelse
// 
//    logprob = logprob - ncp_prior
#if 0
// This should go in a different subclass which corresponds to time-tagged data.
// This version works for the LAT.
  for (deque_t::size_type index = 0; index != num_cells; ++index) {
    double arg = rev_csize[index] - rev_cpop[index] + 1.;
    if (arg > 0.) {
      log_prob[index] = HepStat::gammln(rev_cpop[index] + 1.) + HepStat::gammln(arg) - HepStat::gammln(rev_csize[index] + 2.);
    } else {
      log_prob[index] = eps;
    }
  }
#endif
// 
//  endif else begin                               ;datatype = 3
// 
// ;-----------------------------
// ; posterior for binned data
// ;-----------------------------
// 
//    logprob = lngamma(cell_pops + 1.) - (cell_pops + 1.) * alog(cell_sizes)
// 
//    logprob = logprob - ncp_prior
  for (deque_t::size_type index = 0; index != num_cells; ++index)
    log_prob[index] = HepStat::gammln(rev_cpop[index] + 1.) - (rev_cpop[index] + 1.) * log(rev_csize[index]) - s_ncp_prior;
// 
// ;z = x+1
// ;lnfactx = -z + (z-0.5)*alog(z) + 0.918938533 + 1/(12.*z) - 1/(360.*z^3) + 1/(1260.*z^5) - 1/(1680.*z^7)
// 
//  endelse
////////////////////////////////////////////////////////////////////////////////
}

void TestBurstFitApp::run() {
  using namespace st_facilities;
  m_data_dir = Env::getDataDir("burstFit");

  // Names of test files.
  const char * files[] = { "binData764.fits", "binData1039.fits", "binData1406.fits", "binData2193.fits", "binData2197.fits",
    "binData2387.fits", "binData2665.fits", "binData2711.fits", "binData2863.fits", "binData3256.fits",
    "binData3257.fits", "binData5387.fits", "binData5415.fits", "binData6147.fits", "binData6414.fits",
    "binData6504.fits" };

  // Test them all.
  for (const char ** f_p = files; f_p != files + sizeof(files)/sizeof(const char *); ++f_p)
    testFile(*f_p);
}

void TestBurstFitApp::testFile(const std::string & file_name) {
  using namespace st_facilities;
  using namespace tip;

  // Write file name.
  std::cout << "Test-processing " << file_name << std::endl;

  // Open input file.
  std::auto_ptr<const Table> table(IFileSvc::instance().readTable(Env::appendFileName(m_data_dir, file_name), "1"));

  std::vector<double> cell_pop(table->getNumRecords());
  std::vector<double>::iterator out_itor = cell_pop.begin();
  Table::ConstIterator in_itor = table->begin();
  double cell_size = (*in_itor)["CELLSIZE"].get();
  for (Table::ConstIterator in_itor = table->begin(); in_itor != table->end(); ++in_itor, ++out_itor) {
    *out_itor = (*in_itor)["SUM"].get();
  }

  // Compute blocking needed.
  BayesianBlocker bb(cell_size, cell_pop);

  bb.getBlocks();
}

st_app::StAppFactory<TestBurstFitApp> g_factory;
