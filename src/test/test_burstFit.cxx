/** \file test_burstFit.cxx
    \brief Test application for burstFit.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "burstFit/BurstModel.h"
#include "burstFit/NegativeStat.h"

#include "evtbin/BayesianBinner.h"
#include "evtbin/Hist1D.h"

#include "optimizers/ChiSq.h"
#include "optimizers/Minuit.h"
#include "optimizers/dArg.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_graph/Engine.h"
#include "st_graph/IFrame.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

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
  using namespace burstFit;
  using namespace st_facilities;
  using namespace tip;

  typedef std::vector<double> vec_t;

  // Write file name.
  std::cout << "Test-processing " << file_name << std::endl;

  // Open input file.
  std::auto_ptr<const Table> table(IFileSvc::instance().readTable(Env::appendFileName(m_data_dir, file_name), "1"));

  // Cell populations (original data Y axis).
  vec_t cell_pop(table->getNumRecords());

  // Time domain (original data X axis).
  vec_t domain(cell_pop.size());

  // Input to the Bayesian binner is the original data X values expressed as a sequence of intervals.
  evtbin::BayesianBinner::IntervalCont_t intervals(cell_pop.size());

  // Output index used to populate the several data arrays.
  vec_t::size_type out_index = 0;

  // Iterate over input table.
  double time = 0.;
  for (Table::ConstIterator in_itor = table->begin(); in_itor != table->end(); ++in_itor, ++out_index) {
    domain[out_index] = time;
    time += (*in_itor)["CELLSIZE"].get();
    intervals[out_index] = evtbin::Binner::Interval(domain[out_index], time);
    cell_pop[out_index] = (*in_itor)["SUM"].get();
  }

  // Compute blocking needed.
  evtbin::BayesianBinner bb(intervals, cell_pop.begin(), "SUM");

  // This is based on blocker.pro:
  // tt_blocks = CPs*binscale                  ;convert the CPs to seconds
  // num_blocks = n_elements(tt_blocks) - 1
  // nblocks(itrig) = num_blocks
  // yy_blocks = fltarr(num_blocks)
  // pltfact = pltscale/binscale

  // for iblock = 0,num_blocks-1 do begin      ;intensities in counts/plotting bin
  // yy_blocks(iblock) = total( data(CPs(iblock):CPs(iblock+1)-1) )
  // yy_blocks(iblock) = pltfact * yy_blocks(iblock) / (CPs(iblock+1) - CPs(iblock))
  // endfor
  // yy_blocks = [yy_blocks, 0.]

  // Bin profile into a histogram.
  evtbin::Hist1D hist(bb);
  for (vec_t::size_type index = 0; index != cell_pop.size(); ++index) {
    long block_index = bb.computeIndex(domain[index]);
    
    // Weight for histogram is the cell population multiplied by the original bin width and divided by the bayesian block width.
    hist.fillBin(domain[index],
      cell_pop[index] * (intervals[index].end() - intervals[index].begin()) / bb.getInterval(block_index).width());
  }

  // Create arrays to hold the block domain for plotting purposes.
  vec_t time_start(bb.getNumBins(), 0.);
  vec_t time_stop(bb.getNumBins(), 0.);

  // Loop over the blocks, storing the start/stop time for each.
  for (long block_index = 0; block_index != bb.getNumBins(); ++block_index) {
    evtbin::Binner::Interval interval = bb.getInterval(block_index);

    // Start/stop time is the start/stop point of the interval, scaled by the cell size.
    time_start[block_index] = interval.begin();
    time_stop[block_index] = interval.end();
  }

  // Create a model for this data set.
  BurstModel model(&hist);

  // Compute first peak guess as a function of the domain.
  vec_t guess(domain.size());

  for (vec_t::size_type index = 0; index != domain.size(); ++index) {
    optimizers::dArg arg(domain[index]);
    guess[index] = model.value(arg);
  }

  optimizers::ChiSq chi_sq(domain, cell_pop, &model);
  NegativeStat stat(&chi_sq);

  optimizers::Minuit opt(stat);

  std::vector<double> coeff;
  model.getParamValues(coeff);
  std::clog << "Before fit:" << std::endl;
  for (std::vector<double>::iterator itor = coeff.begin(); itor != coeff.end(); ++itor) {
    std::clog << "\tPar " << itor - coeff.begin() << " is " << *itor << std::endl;
  }
  std::clog << "\tChiSq is " << chi_sq.value() / chi_sq.dof() << std::endl;

  bool converged = false;
  try {
    opt.find_min();
    converged = true;
  } catch (const std::exception & x) {
    std::clog << x.what() << std::endl;
  }

  std::clog << "After fit:" << std::endl;
  coeff.clear();
  model.getParamValues(coeff);
  for (std::vector<double>::iterator itor = coeff.begin(); itor != coeff.end(); ++itor) {
    std::clog << "\tPar " << itor - coeff.begin() << " is " << *itor << std::endl;
  }
  std::clog << "\tChiSq is " << chi_sq.value() / chi_sq.dof() << std::endl;

  vec_t fit(domain.size());
  for (vec_t::size_type index = 0; index != domain.size(); ++index) {
    optimizers::dArg arg(domain[index]);
    fit[index] = model.value(arg);
  }

  // Plot blocks and data:
  try {
    // Get plot engine.
    st_graph::Engine & engine(st_graph::Engine::instance());

    // Create main frame.
    std::auto_ptr<st_graph::IFrame> mf(engine.createMainFrame(0, 600, 400, "test_burstFit"));

    // Create plot frame.
    std::auto_ptr<st_graph::IFrame> pf(engine.createPlotFrame(mf.get(), file_name, 600, 400));

    // Plot the data.
    std::auto_ptr<st_graph::IPlot> data_plot(engine.createPlot(pf.get(), "hist",
      st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
      st_graph::PointSequence<vec_t::iterator>(cell_pop.begin(), cell_pop.end())));

    // Overplot the blocks.
#if 0
    std::auto_ptr<st_graph::IPlot> block_plot(engine.createPlot(pf.get(), "hist",
      st_graph::IntervalSequence<vec_t::iterator>(time_start.begin(), time_start.end(), time_stop.begin()),
      st_graph::PointSequence<evtbin::Hist1D::ConstIterator>(hist.begin(), hist.end())));
#endif

#if 0
    // Overplot the initial guess.
    std::auto_ptr<st_graph::IPlot> guess_plot(engine.createPlot(pf.get(), "hist",
      st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
      st_graph::PointSequence<vec_t::iterator>(guess.begin(), guess.end())));
#endif

    // Overplot the final fit if it converged.
//    if (converged) {
      std::auto_ptr<st_graph::IPlot> fit_plot(engine.createPlot(pf.get(), "hist",
        st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
        st_graph::PointSequence<vec_t::iterator>(fit.begin(), fit.end())));
//    }

    engine.run();
  } catch (const std::exception & x) {
    std::clog << "Could not display test plots:" << x.what() << std::endl;
    // Ignore errors with the plotting.
  }

}

st_app::StAppFactory<TestBurstFitApp> g_factory;
