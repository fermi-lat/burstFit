/**
    \mainpage burstFit package

    \author  James Peachey peachey@milkyway.gsfc.nasa.gov
             Lawrence Brown elwin@milkyway.gsfc.nasa.gov

    \section synopsis Synopsis
    The burstFit package contains a library and an application,
    gtburstfit. The application analyzes burst light curves
    by applying a Bayesian algorithm to determine the optimum
    set of blocks to follow the burst profile shape. These
    blocks may also be used by gtburstfit to construct a
    pulse model, which gtburstfit can fit to the data.

    \section parameters Application Parameters

    \subsection key Key To Parameter Descriptions
\verbatim
Automatic parameters:
par_name [ = value ] type

Hidden parameters:
(par_name = value ) type

Where "par_name" is the name of the parameter, "value" is the
default value, and "type" is the type of the parameter. The
type is enclosed in square brackets.

Examples:
infile [file]
    Describes an automatic (queried) file-type parameter with
    no default value.

(plot = yes) [bool]
    Describes a hidden bool-type parameter named plot, whose
    default value is yes (true).
\endverbatim

    <a name="gtburstfit_parameters"></a>
    \subsection gtburstfit_parameters gtburstfit Parameters

\verbatim
evfile [file]
    Name of input event file. Currently, only binned events
    (binned light curves in a RATE FITS extension) are
    supported, for both GBM and LAT data.

(evtable = RATE) [string]
    Name of table containing event data.

(fit = YES) [bool]
    If yes, causes gtburstfit to fit parameters to minimize
    Chi Square for the deviation from the fit model. If no,
    only the definitions of Bayesian blocks found by the
    tool will be presented.

(plot = YES) [bool]
    If yes, causes gtburstfit to display a plot showing
    the data, the Bayesian Blocks which were created for
    the data, and the results of the fit. (The fit will
    only be displayed if the fit parameter is set to YES).

(ncpprior = 9.) [double]
    This number represents a heuristic bias against creating
    new blocks, used to allow some control over the number
    of blocks created. If too many (too few) blocks are
    being created for a given set of data, the user may
    increase (decrease) this parameter value.

\endverbatim

    \section algorithm Algorithms
    The Bayesian method used to define the blocks was developed
    by Jeff Scargle and Jay Norris. The model used to fit
    the data is a constant background added to a sum of terms
    of the form:
\verbatim

    A * exp(tau1 / (t - t0) + (t - t0) / tau2)

\endverbatim

    The number of these terms is determined automatically by
    gtburstfit by looking for peaks and valleys in the Bayesian
    block definitions, and using these to determine the number of
    peaks detected, and their approximate positions, heights
    and widths.
    
*/
