/** \file test_burstFit.cxx
    \brief Test application for burstFit.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "burstFit/BurstModel.h"
#include "burstFit/NegativeStat.h"

#include "evtbin/BayesianBinner.h"
#include "evtbin/Hist1D.h"

#include "optimizers/ChiSq.h"
#include "optimizers/Minuit.h"
#include "optimizers/dArg.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_graph/Engine.h"
#include "st_graph/IFrame.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace {
  std::ostream & operator <<(std::ostream & os, const burstFit::BurstModel & model) {
    st_stream::OStream st_stream_os(false);
    st_stream_os.connect(os);
    st_stream_os << model;
    return os;
  }
}

class TestBurstFitApp : public st_app::StApp {
  public:
    virtual ~TestBurstFitApp() throw() {}

    virtual void run();

    virtual void testFile(const std::string & file_name, bool plot);

  private:
    std::string m_data_dir;
};

void TestBurstFitApp::run() {
  using namespace st_facilities;
  m_data_dir = Env::getDataDir("burstFit");

  // Prompt for parameters.
  st_app::AppParGroup & pars(getParGroup());
  pars.Prompt();

  bool plot = pars["plot"];

  // Names of test files.
  const char * files[] = { "ft1tiny.fits", "n0_p05.lc", "binData764.fits", "binData1039.fits",
    "binData1406.fits", "binData2193.fits", "binData2197.fits", "binData2387.fits", "binData2665.fits",
    "binData2711.fits", "binData2863.fits", "binData3256.fits", "binData3257.fits", "binData5387.fits",
    "binData5415.fits", "binData6147.fits", "binData6414.fits", "binData6504.fits" };

  // Test them all.
  for (const char ** f_p = files; f_p != files + sizeof(files)/sizeof(const char *); ++f_p)
    testFile(*f_p, plot);
}

void TestBurstFitApp::testFile(const std::string & file_name, bool plot) {
  using namespace burstFit;
  using namespace st_facilities;
  using namespace tip;

  typedef std::vector<double> vec_t;

  // Write file name.
  std::clog << "Test-processing " << file_name << std::endl;

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

  // Determine whether there is a "TIME" field.
  std::string time_field;
  bool have_time_field = false;
  try {
    table->getFieldIndex("TIME");
    time_field = "TIME";
    have_time_field = true;
  } catch (const std::exception &) {
  }

  // Determine names of fields containing the cell size and content.
  std::string cell_size_field;
  bool have_cell_size_field = false;
  try {
    table->getFieldIndex("TIMEDEL");
    cell_size_field = "TIMEDEL";
    have_cell_size_field = true;
  } catch (const std::exception &) {
  }
  if (!have_cell_size_field) {
    try {
      table->getFieldIndex("CELLSIZE");
      cell_size_field = "CELLSIZE";
      have_cell_size_field = true;
    } catch (const std::exception &) {
    }
  }

  if (!have_time_field && !have_cell_size_field)
    throw std::runtime_error("Could not find field TIME, TIMEDEL or CELLSIZE in file " + file_name);

  std::string cell_pop_field;
  bool have_cell_pop_field = false;
  try {
    table->getFieldIndex("COUNTS");
    cell_pop_field = "COUNTS";
    have_cell_pop_field = true;
  } catch (const std::exception &) {
  }
  if (!have_cell_pop_field) {
    try {
      table->getFieldIndex("SUM");
      cell_pop_field = "SUM";
      have_cell_pop_field = true;
    } catch (const std::exception &) {
    }
  }

  // Iterate over input table, extract data.
  double time = 0.;
  if (have_time_field) time = (*table->begin())[time_field].get();

  for (Table::ConstIterator in_itor = table->begin(); in_itor != table->end(); ++in_itor, ++out_index) {
    Table::ConstIterator next_itor = in_itor;
    ++next_itor;
    domain[out_index] = time;
    if (have_time_field && (next_itor != table->end())) {
      time = (*(next_itor))[time_field].get(); 
    } else if (have_cell_size_field) {
      time += (*in_itor)[cell_size_field].get();
    }
    intervals[out_index] = evtbin::Binner::Interval(domain[out_index], time);
    if (have_cell_pop_field) {
      cell_pop[out_index] = (*in_itor)[cell_pop_field].get();
    } else {
      cell_pop[out_index] = 1.;
    }
  }

  // If cell sizes were not explicitly supplied, the last data point must be discarded.
  if (!have_cell_size_field) {
    std::vector<double>::size_type new_size = cell_pop.size() - 1;
    domain.resize(new_size);
    cell_pop.resize(new_size);
    intervals.resize(new_size);
  }

  // Compute blocking needed.
  evtbin::BayesianBinner bb(intervals, cell_pop.begin(), cell_pop_field);

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

  // Create optimizing objective function.
  optimizers::ChiSq chi_sq(domain, cell_pop, &model);

  // The function to MAXIMIZE is the negative of the chi_sq.
  NegativeStat stat(&chi_sq);

  // Create optimizer for the objective function.
  optimizers::Minuit opt(stat);

  // Display fit parameters before fit.
  std::vector<double> coeff;
  model.getParamValues(coeff);
  std::clog << "After initial guess but before any fitting, reduced chi square is " << chi_sq.value() / chi_sq.dof() << std::endl;
  std::clog << "Parameters are:" << std::endl;
  std::clog << model << "\n" << std::endl;

  bool converged = false;
  try {
    opt.find_min();
    converged = true;
  } catch (const std::exception & x) {
    std::clog << x.what() << std::endl;
  }

  coeff.clear(); // Just in case things are malfunctioning badly.
  model.getParamValues(coeff);

  std::clog << "After fit, reduced chi square is " << chi_sq.value() / chi_sq.dof() << std::endl;
  std::clog << "Parameters are:" << std::endl;
  std::clog << model << "\n" << std::endl;

  vec_t fit(domain.size());
  for (vec_t::size_type index = 0; index != domain.size(); ++index) {
    optimizers::dArg arg(domain[index]);
    fit[index] = model.value(arg);
  }

  if (plot) {
    // Plot blocks and/or data:
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
      std::auto_ptr<st_graph::IPlot> block_plot(engine.createPlot(pf.get(), "hist",
        st_graph::IntervalSequence<vec_t::iterator>(time_start.begin(), time_start.end(), time_stop.begin()),
        st_graph::PointSequence<evtbin::Hist1D::ConstIterator>(hist.begin(), hist.end())));
  
      // Overplot the final fit.
      std::auto_ptr<st_graph::IPlot> fit_plot(0);
      if (converged)
        fit_plot.reset(engine.createPlot(pf.get(), "hist",
          st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
          st_graph::PointSequence<vec_t::iterator>(fit.begin(), fit.end())));
  
      engine.run();
    } catch (const std::exception & x) {
      std::clog << "Could not display test plots:" << x.what() << std::endl;
      // Ignore errors with the plotting.
    }
  }

}

st_app::StAppFactory<TestBurstFitApp> g_factory("test_burstFit");
