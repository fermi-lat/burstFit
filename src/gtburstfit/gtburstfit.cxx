/** \file gtburstfit.cxx
    \brief Main burstFit application.
    \author James Peachey, HEASARC/GSSC
*/
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

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

class BurstFitApp : public st_app::StApp {
  public:
    BurstFitApp();
    virtual ~BurstFitApp() throw() {}

    virtual void run();

  private:
    st_stream::StreamFormatter m_os;
};

BurstFitApp::BurstFitApp(): m_os(getName(), "BurstFitApp", 2) {
}

void BurstFitApp::run() {
  // Prompt for parameters.
  st_app::AppParGroup & pars(getParGroup());
  pars.Prompt();
  pars.Save();

  std::string ev_file = pars["evfile"];
  std::string ev_table = pars["evtable"];
  bool fit = pars["fit"];
  bool plot = pars["plot"];
  double ncp_prior = pars["ncpprior"];

  using namespace burstFit;
  using namespace tip;

  typedef std::vector<double> vec_t;

  // Open input file.
  std::auto_ptr<const Table> table(IFileSvc::instance().readTable(ev_file, ev_table));

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
    throw std::runtime_error("Could not find field TIME, TIMEDEL or CELLSIZE in file " + ev_file);

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

  // Compute Bayesian blocks
  evtbin::BayesianBinner bb(intervals, cell_pop.begin(), cell_pop_field, ncp_prior);

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

  vec_t fit_result(domain.size());
  if (fit) {
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
    m_os.info() << "After initial guess but before any fitting, reduced chi square is " << chi_sq.value() / chi_sq.dof() <<
      std::endl;
    m_os.info() << "Parameters are:" << std::endl;
    m_os.info() << model << std::endl << std::endl;
  
    bool converged = false;
    try {
      opt.find_min();
      converged = true;
    } catch (const std::exception & x) {
      m_os.err() << x.what() << std::endl;
    }
  
    coeff.clear(); // Just in case things are malfunctioning badly.
    model.getParamValues(coeff);
  
    m_os.info() << "After fit, reduced chi square is " << chi_sq.value() / chi_sq.dof() << std::endl;
    m_os.info() << "Parameters are:" << std::endl;
    m_os.info() << model << std::endl << std::endl;
  
    for (vec_t::size_type index = 0; index != domain.size(); ++index) {
      optimizers::dArg arg(domain[index]);
      fit_result[index] = model.value(arg);
    }
  } else {
    m_os.info() << "Bayesian Blocks computed for this data set are:" << std::endl;
    m_os.info() << "Interval           Average Counts" << std::endl;
    for (long index = 0; index != bb.getNumBins(); ++index) {
      m_os.info() << "[" << time_start[index] << ", " << time_stop[index] << "]    " << hist[index] << std::endl;
    }
  }

  if (plot) {
    // Plot blocks and/or data:
    try {
      // Get plot engine.
      st_graph::Engine & engine(st_graph::Engine::instance());
  
      // Create main frame.
      std::auto_ptr<st_graph::IFrame> mf(engine.createMainFrame(0, 600, 400, "gtburstfit"));
  
      // Create plot frame.
      std::string::size_type begin = ev_file.find_last_of("/\\");
      if (std::string::npos == begin) begin = 0; else ++begin;
      std::string title = ev_file.substr(begin) + ": Black - data, Red - Bayesian Blocks" + (fit ? ", Green - fit model" : "");
      std::auto_ptr<st_graph::IFrame> pf(engine.createPlotFrame(mf.get(), title, 600, 400));
  
      // Plot the data.
      std::auto_ptr<st_graph::IPlot> data_plot(engine.createPlot(pf.get(), "hist",
        st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
        st_graph::PointSequence<vec_t::iterator>(cell_pop.begin(), cell_pop.end())));
  
      // Overplot the blocks.
      std::auto_ptr<st_graph::IPlot> block_plot(engine.createPlot(pf.get(), "hist",
        st_graph::IntervalSequence<vec_t::iterator>(time_start.begin(), time_start.end(), time_stop.begin()),
        st_graph::PointSequence<evtbin::Hist1D::ConstIterator>(hist.begin(), hist.end())));
  
      std::auto_ptr<st_graph::IPlot> fit_plot(0);
      if (fit) {
        // Overplot the final fit.
        fit_plot.reset(engine.createPlot(pf.get(), "hist",
          st_graph::LowerBoundSequence<vec_t::iterator>(domain.begin(), domain.end()),
          st_graph::PointSequence<vec_t::iterator>(fit_result.begin(), fit_result.end())));
      }
  
      engine.run();
    } catch (const std::exception & x) {
      m_os.err() << "Could not display plot:" << x.what() << std::endl;
      // Ignore errors with the plotting.
    }
  }

}

st_app::StAppFactory<BurstFitApp> g_factory("gtburstfit");
