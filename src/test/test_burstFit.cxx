/** \file test_burstFit.cxx
    \brief Test application for burstFit.
    \authors Lawrence Brown, HEASARC
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "burstFit/BayesianBlocker.h"

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
  burstFit::BayesianBlocker bb(cell_size, cell_pop);

  bb.getBlocks();
}

st_app::StAppFactory<TestBurstFitApp> g_factory;
