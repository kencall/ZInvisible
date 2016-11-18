IDIR       = .
ODIR       = obj
SDIR       = .
RSDIR      = $(CMSSW_BASE)/src/SusyAnaTools/Tools

TTIDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/include
TPIDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/include
TPTDIR     = $(CMSSW_BASE)/src/TopTagger/Tools

CXX        = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -std=c++0x -I$(LHAPDF_DATA_PATH)/../../include
## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

LIBS       = $(shell root-config --glibs)
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -lTopTaggerTopTagger -lTopTaggerCfgParser
LHAPDFLIB  = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF

#OBJS       = $(patsubst %, $(ODIR)/%, $(OBJ))

PROGRAMS = calcEff harvestHists makePlots systematics nJetWgtSyst moneyplot plot2d scalefactorplot systplot miniMakePlots beffCalc makeSignalHistograms makeTopPlots

all: mkobj PyWraps $(PROGRAMS)


mkobj:
	@mkdir -p obj

PyWraps: $(ODIR)/samplesModule.so $(ODIR)/ScaleFactors.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(RSDIR)/samples.cc $(RSDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(RSDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/ScaleFactors.so: $(ODIR)/ScaleFactorsPyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/ScaleFactorsPyWrap.o: $(SDIR)/ScaleFactors.cc $(SDIR)/ScaleFactors.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(TPTDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(TPTDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

calcEff: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/calcEff.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/Plotter.o $(ODIR)/MiniTupleMaker.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

harvestHists: $(ODIR)/harvestHists.o
	$(LD) $^ $(LIBS) -o $@

makePlots: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o $(ODIR)/MiniTupleMaker.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

makeTopPlots: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/Plotter.o $(ODIR)/MakeTopPlots.o $(ODIR)/MiniTupleMaker.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

miniMakePlots: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/miniTupleSamples.o $(ODIR)/Plotter.o $(ODIR)/MinitupleMakePlots.o $(ODIR)/MiniTupleMaker.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

systematics: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/miniTupleSamples.o $(ODIR)/Systematic.o  $(ODIR)/Plotter.o $(ODIR)/systematics.o $(ODIR)/MiniTupleMaker.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

nJetWgtSyst: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/nJetWgtSyst.o $(ODIR)/miniTupleSamples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

moneyplot: $(ODIR)/moneyplot.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -o $@

scalefactorplot: $(ODIR)/scalefactorplot.o
	$(LD) $^ $(LIBS) -o $@

systplot: $(ODIR)/systplot.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -o $@

plot2d: $(ODIR)/plot2D.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/Plotter.o $(ODIR)/PDFUncertainty.o $(ODIR)/RegisterFunctions.o $(ODIR)/Systematic.o $(ODIR)/samples.o $(ODIR)/MiniTupleMaker.o $(ODIR)/searchBins.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/PlotUtility.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(LHAPDFLIB) -o $@

beffCalc: $(ODIR)/bTagEfficiencyCalc.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) -o $@

makeSignalHistograms: $(ODIR)/makeSignalHistograms.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o 
	$(LD) $^ $(LIBS) $(MT2LIB) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
