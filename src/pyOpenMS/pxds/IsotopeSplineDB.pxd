from Types cimport *
from libcpp cimport bool
from String cimport *
from Map cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/IsotopeSplineDB.h>" namespace "OpenMS":

    cdef cppclass IsotopeSplineDB "OpenMS::IsotopeSplineDB":
        # wrap-manual-memory

        IsotopeSplineDB(IsotopeSplineDB) nogil except + #wrap-ignore

        # Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        IsotopeDistribution estimateFromPeptideWeight(double average_weight, UInt max_depth) nogil except +

        # Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs
        IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S, UInt max_depth) nogil except +

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
        # and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +


        bool inModelBounds(double average_weight, UInt max_depth) nogil except +

        bool inModelBounds(double average_weight, UInt S, UInt max_depth) nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/IsotopeSplineDB.h>" namespace "OpenMS::IsotopeSplineDB":

    const IsotopeSplineDB* getInstance() nogil except + # wrap-ignore

