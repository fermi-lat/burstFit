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

#include "burstFit/BayesianBinner.h"
#include "burstFit/BayesianBlocker.h"

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
  using namespace st_facilities;
  using namespace tip;

  // Write file name.
  std::cout << "Test-processing " << file_name << std::endl;

  // Open input file.
  std::auto_ptr<const Table> table(IFileSvc::instance().readTable(Env::appendFileName(m_data_dir, file_name), "1"));

  std::vector<double> cell_pop(table->getNumRecords());
  std::vector<double>::iterator out_itor = cell_pop.begin();
  std::vector<double> domain(cell_pop.size());
  std::vector<double>::iterator domain_itor = domain.begin();
  Table::ConstIterator in_itor = table->begin();
  double cell_size = (*in_itor)["CELLSIZE"].get();
  for (Table::ConstIterator in_itor = table->begin(); in_itor != table->end(); ++in_itor, ++out_itor, ++domain_itor) {
    *domain_itor = (out_itor - cell_pop.begin()) * cell_size;
    *out_itor = (*in_itor)["SUM"].get();
  }

  // Compute blocking needed.
  evtbin::BayesianBinner bb("SUM", cell_size, cell_pop.begin(), cell_pop.end());

  typedef std::vector<double> vec_t;

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

  // Fractional threshold used to determine significant deviations from background.
  double threshold_fract = .05;

  // Create arrays to hold the blocks.
  // These play the role of tt_blocks and yy_blocks in pulsefitter6.pro.
  vec_t time_start(bb.getNumBins(), 0.);
  vec_t time_stop(bb.getNumBins(), 0.);
  vec_t block(bb.getNumBins(), 0.);

  // For each cell, bin it into a block, using the width of the block to compute the mean.
  for (vec_t::size_type index = 0; index != cell_pop.size(); ++index) {
    // Determine to which block this cell belongs.
    long block_index = bb.computeIndex(index);

    // Add the cell population to this block. Note: block holds the total number of counts at this point.
    block[block_index] += cell_pop[index];
  }
  
  // Loop over the blocks, computing the time for each block and also computing mean for each block from total counts.
  for (long block_index = 0; block_index != bb.getNumBins(); ++block_index) {
    evtbin::Binner::Interval interval = bb.getInterval(block_index);

    // Start/stop time is the start/stop point of the interval, scaled by the cell size.
    time_start[block_index] = cell_size * interval.begin();
    time_stop[block_index] = cell_size * interval.end();

    // Convert block array into the mean number of counts by dividing by the block width.
    block[block_index] /= interval.width();
  }

  // The following is based on findpvs3.pro:
  if (block.size() >= 3) {
    // First and last block are used to set a background level, used to determine thresholds.
    double background = (block.front() + block.back()) / 2.;

    // Also get maximum block value, which plays a role in determining whether threshold is exceeded.
    double max_block = *std::max_element(block.begin(), block.end());

    // TODO: Double check this with experts: this test seems like it will always pass. Compare to findpvs3.pro.
    // Threshold for a PEAK counting as "significant" is that it exceed background by > threshold based on
    // the maximum BLOCK value. (Note: maximum peak != maximum block in general.)
    double threshold = threshold_fract * (max_block - background);

    // First look for peaks, which are by definition any/all local maxima.
    // Container to hold the indices of peak locations.
    std::vector<vec_t::size_type> peak_index;
    bool significant_peaks = false;
    double total_peak = 0.;
    for (vec_t::size_type index = 1; index != block.size() - 1; ++index) {
      if (block[index - 1] < block[index] && block[index] > block[index + 1]) {
        // All local maxima count as peaks.
        peak_index.push_back(index);

        // Only some maxima count as significant for purposes of determining whether to continue.
        if (block[index] - background > threshold) significant_peaks = true;

        // Compute the sum of all peaks for use below.
        total_peak += block[index];
      }
    }

    if (significant_peaks) {
      // Container to hold the indices of valley locations.
      std::vector<vec_t::size_type> valley_index;

      // Determine the average peak height for use in setting the threshold for valleys.
      double average_peak = total_peak / peak_index.size() - background;

      // Threshold is a fraction of the average peak height.
      threshold = threshold_fract * average_peak;

      // Special handling for the 0th valley.
      // Look for first rise; skip all blocks before the first peak which are below threshold.
      vec_t::size_type index;
      for (index = 1; index != peak_index.front() && block[index] - block.front() <= threshold; ++index) {}

      // 0th valley is the block before the first rise.
      valley_index.push_back(index - 1);

      // Handle interjacent valleys, which are minima between the maxima.
      for (std::vector<vec_t::size_type>::iterator itor = peak_index.begin(); itor != peak_index.end() - 1; ++itor) {
        valley_index.push_back(std::min_element(block.begin() + *itor, block.begin() + *(itor + 1)) - block.begin());
      }

      // Special handling for the last valley.
      // Look for the last decay, that is the first block after the last peak which is below threshold.
      for (index = peak_index.back() + 1; index != block.size() - 1 && block[index] - block.back() > threshold; ++index) {}

      // Last valley is the block after the last decay.
      valley_index.push_back(index);
    } else {
std::clog << "No significant peaks." << std::endl;
    }

  } else {
std::clog << "Fewer than 3 blocks." << std::endl;
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
    std::auto_ptr<st_graph::IPlot> block_plot(engine.createPlot(pf.get(), "hist",
      st_graph::IntervalSequence<vec_t::iterator>(time_start.begin(), time_start.end(), time_stop.begin()),
      st_graph::PointSequence<vec_t::iterator>(block.begin(), block.end())));

    engine.run();
  } catch (const std::exception &) {
    std::clog << "Could not display plots." << std::endl;
    // Ignore errors with the plotting.
  }

}

st_app::StAppFactory<TestBurstFitApp> g_factory;
